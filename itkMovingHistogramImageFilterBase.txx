/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMovingHistogramImageFilterBase.txx,v $
  Language:  C++
  Date:      $Date: 2004/04/30 21:02:03 $
  Version:   $Revision: 1.14 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMovingHistogramImageFilterBase_txx
#define __itkMovingHistogramImageFilterBase_txx

#include "itkMovingHistogramImageFilterBase.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkOffset.h"
#include "itkProgressReporter.h"
#include "itkNumericTraits.h"

#ifndef zigzag

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageLinearConstIteratorWithIndex.h"

#endif

namespace itk {


template<class TInputImage, class TOutputImage, class TKernel>
MovingHistogramImageFilterBase<TInputImage, TOutputImage, TKernel>
::MovingHistogramImageFilterBase()
  : m_Kernel()
{
  m_PixelsPerTranslation = 0;
  this->SetRadius( 1 );
}


template<class TInputImage, class TOutputImage, class TKernel>
void
MovingHistogramImageFilterBase<TInputImage, TOutputImage, TKernel>
::SetRadius( const RadiusType & radius )
{
  KernelType kernel;
  kernel.SetRadius( radius );
  for( typename KernelType::Iterator kit=kernel.Begin(); kit!=kernel.End(); kit++ )
    {
    *kit = 1;
    }
  this->SetKernel( kernel );
}


template<class TInputImage, class TOutputImage, class TKernel>
void
MovingHistogramImageFilterBase<TInputImage, TOutputImage, TKernel>
::SetRadius( unsigned long radius )
{
  RadiusType rad;
  rad.Fill( radius );
  this->SetRadius( rad );
}


template<class TInputImage, class TOutputImage, class TKernel>
void
MovingHistogramImageFilterBase<TInputImage, TOutputImage, TKernel>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  
  // get pointers to the input and output
  typename Superclass::InputImagePointer  inputPtr = 
    const_cast< TInputImage * >( this->GetInput() );
  
  if ( !inputPtr )
    {
    return;
    }

  // get a copy of the input requested region (should equal the output
  // requested region)
  typename TInputImage::RegionType inputRequestedRegion;
  inputRequestedRegion = inputPtr->GetRequestedRegion();

  // pad the input requested region by the operator radius
  inputRequestedRegion.PadByRadius( m_Kernel.GetRadius() );

  // crop the input requested region at the input's largest possible region
  if ( inputRequestedRegion.Crop(inputPtr->GetLargestPossibleRegion()) )
    {
    inputPtr->SetRequestedRegion( inputRequestedRegion );
    return;
    }
  else
    {
    // Couldn't crop the region (requested region is outside the largest
    // possible region).  Throw an exception.

    // store what we tried to request (prior to trying to crop)
    inputPtr->SetRequestedRegion( inputRequestedRegion );
    
    // build an exception
    InvalidRequestedRegionError e(__FILE__, __LINE__);
    OStringStream msg;
    msg << static_cast<const char *>(this->GetNameOfClass())
        << "::GenerateInputRequestedRegion()";
    e.SetLocation(msg.str().c_str());
    e.SetDescription("Requested region is (at least partially) outside the largest possible region.");
    e.SetDataObject(inputPtr);
    throw e;
    }
}

template<class TInputImage, class TOutputImage, class TKernel>
void
MovingHistogramImageFilterBase<TInputImage, TOutputImage, TKernel>
::SetKernel( const KernelType& kernel )
{
  // first, build the list of offsets of added and removed pixels when the 
  // structuring element move of 1 pixel on 1 axis; do it for the 2 directions
  // on each axes.
  
  // transform the structuring element in an image for an easier
  // access to the data
  typedef Image< bool, TInputImage::ImageDimension > BoolImageType;
  typename BoolImageType::Pointer tmpSEImage = BoolImageType::New();
  tmpSEImage->SetRegions( kernel.GetSize() );
  tmpSEImage->Allocate();
  RegionType tmpSEImageRegion = tmpSEImage->GetRequestedRegion();
  ImageRegionIteratorWithIndex<BoolImageType> kernelImageIt;
  kernelImageIt = ImageRegionIteratorWithIndex<BoolImageType>(tmpSEImage, tmpSEImageRegion);
  kernelImageIt.GoToBegin();
  KernelIteratorType kernel_it = kernel.Begin();
  OffsetListType kernelOffsets;

  // create a center index to compute the offset
  IndexType centerIndex;
  for( unsigned axis=0; axis<ImageDimension; axis++)
    { centerIndex[axis] = kernel.GetSize()[axis] / 2; }
  
  unsigned long count = 0;
  while( !kernelImageIt.IsAtEnd() )
    {
    kernelImageIt.Set( *kernel_it > 0 );
    if( *kernel_it > 0 )
      {
      kernelImageIt.Set( true );
      kernelOffsets.push_front( kernelImageIt.GetIndex() - centerIndex );
      count++;
      }
    else
      { kernelImageIt.Set( false ); }
    ++kernelImageIt;
    ++kernel_it;
    }
    

  // verify that the kernel contain at least one point
  if( count == 0 )
    { itkExceptionMacro( << "The kernel must contain at least one point." ); }

  // no attribute should be modified before here to avoid setting the filter in a bad status
  // store the kernel !!
  m_Kernel = kernel;

  // clear the already stored values
  m_AddedOffsets.clear();
  m_RemovedOffsets.clear();
  //m_Axes

  // store the kernel offset list
  m_KernelOffsets = kernelOffsets;

  typename itk::FixedArray< unsigned long, ImageDimension > axisCount;
  axisCount.Fill( 0 );

  for( unsigned axis=0; axis<ImageDimension; axis++)
    {
    OffsetType refOffset;
    refOffset.Fill( 0 );
    for( int direction=-1; direction<=1; direction +=2)
      {
      refOffset[axis] = direction;
      for( kernelImageIt.GoToBegin(); !kernelImageIt.IsAtEnd(); ++kernelImageIt)
        {
        IndexType idx = kernelImageIt.GetIndex();
        
        if( kernelImageIt.Get() )
          {
          // search for added pixel during a translation
          IndexType nextIdx = idx + refOffset;
          if( tmpSEImageRegion.IsInside( nextIdx ) )
            {
            if( !tmpSEImage->GetPixel( nextIdx ) )
              {
                m_AddedOffsets[refOffset].push_front( nextIdx - centerIndex );
                axisCount[axis]++;
              }
            }
          else
            {
              m_AddedOffsets[refOffset].push_front( nextIdx - centerIndex );
              axisCount[axis]++;
            }
          // search for removed pixel during a translation
          IndexType prevIdx = idx - refOffset;
          if( tmpSEImageRegion.IsInside( prevIdx ) )
            {
            if( !tmpSEImage->GetPixel( prevIdx ) )
              {
                m_RemovedOffsets[refOffset].push_front( idx - centerIndex );
                axisCount[axis]++;
              }
            }
          else
            {
              m_RemovedOffsets[refOffset].push_front( idx - centerIndex );
              axisCount[axis]++;
            }

          }
        }
      }
    }
    
    // search for the best axis
    typedef typename std::set<DirectionCost> MapCountType;
    MapCountType invertedCount;
    for( unsigned i=0; i<ImageDimension; i++ )
      {
      invertedCount.insert( DirectionCost( i, axisCount[i] ) );
      }

    int i=0;
    for( typename MapCountType::iterator it=invertedCount.begin(); it!=invertedCount.end(); it++, i++)
      {
      m_Axes[i] = it->m_Dimension;
      }

    m_PixelsPerTranslation = axisCount[m_Axes[ImageDimension - 1]] / 2;  // divided by 2 because there is 2 directions on the axis
}


template<class TInputImage, class TOutputImage, class TKernel>
void
MovingHistogramImageFilterBase<TInputImage, TOutputImage, TKernel>
::GetDirAndOffset(const IndexType LineStart, 
		  const IndexType PrevLineStart,
		  const int ImageDimension,
		  OffsetType &LineOffset,
		  OffsetType &Changes,
		  int &LineDirection)
{
  // when moving between lines in the same plane there should be only
  // 1 non zero (positive) entry in LineOffset.
  // When moving between planes there will be some negative ones too.
  LineOffset = Changes = LineStart - PrevLineStart;
  for (int y=0;y<ImageDimension;y++) 
    {
    if (LineOffset[y] > 0)
      {
      LineOffset[y]=1;  // should be 1 anyway
      LineDirection=y;
      }
    else
      {
      LineOffset[y]=0;
      }
    }
}


template<class TInputImage, class TOutputImage, class TKernel>
void
MovingHistogramImageFilterBase<TInputImage, TOutputImage, TKernel>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Kernel: " << m_Kernel << std::endl;
}

}// end namespace itk
#endif
