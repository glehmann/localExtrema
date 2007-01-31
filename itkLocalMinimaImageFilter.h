/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLocalMinimaImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2004/04/30 21:02:03 $
  Version:   $Revision: 1.15 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkLocalMinimaImageFilter_h
#define __itkLocalMinimaImageFilter_h

#include "itkMovingHistogramImageFilter.h"
#include <map>

namespace itk {

namespace Function {
template <class TInputPixel, class TOutputPixel>
class LocalMinimaHistogram
{
public:
  LocalMinimaHistogram()
    {
    if( useVectorBasedAlgorithm() )
      { initVector(); }
    }
  ~LocalMinimaHistogram(){}

  LocalMinimaHistogram * Clone()
    {
    LocalMinimaHistogram * result = new LocalMinimaHistogram();
    result->m_Map = this->m_Map;
    result->m_Vector = this->m_Vector;
    result->m_Min = this->m_Min;
    result->m_Max = this->m_Max;
    result->m_Count = this->m_Count;
    result->m_ForegroundValue = this->m_ForegroundValue;
    result->m_BackgroundValue = this->m_BackgroundValue;
    result->m_FlatToForeground = this->m_FlatToForeground;
    result->m_Threshold = this->m_Threshold;
    return result;
    }

  inline void AddBoundary() {}

  inline void RemoveBoundary() {}

  inline void AddPixel( const TInputPixel &p )
    {
    if( useVectorBasedAlgorithm() )
      { AddPixelVector( p ); }
    else
      { AddPixelMap( p ); }
    }

  inline void RemovePixel( const TInputPixel &p )
    {
    if( useVectorBasedAlgorithm() )
      { RemovePixelVector( p ); }
    else
      { RemovePixelMap( p ); }
    }

  inline const TOutputPixel & GetValue( const TInputPixel & centerPixel )
    {
    if( useVectorBasedAlgorithm() )
      { return GetValueVector( centerPixel ); }
    else
      { return GetValueMap( centerPixel ); }
    }


  static inline bool useVectorBasedAlgorithm()
    {
    // bool, short and char are acceptable for vector based algorithm: they do not require
    // too much memory. Other types are not usable with that algorithm
    return typeid(TInputPixel) == typeid(unsigned char)
        || typeid(TInputPixel) == typeid(signed char)
        || typeid(TInputPixel) == typeid(unsigned short)
        || typeid(TInputPixel) == typeid(signed short)
        || typeid(TInputPixel) == typeid(bool);
    }







  //
  // the map based algorithm
  //

  typedef typename std::map< TInputPixel, unsigned long > MapType;

  inline void AddPixelMap( const TInputPixel &p )
    { m_Map[ p ]++; }

  inline void RemovePixelMap( const TInputPixel &p )
    { m_Map[ p ]--; }

  inline const TOutputPixel & GetValueMap( const TInputPixel & centerPixel )
    {
    if( centerPixel > m_Threshold )
      {
      return m_BackgroundValue;
      }

    // clean the map
    typename MapType::iterator mapIt = m_Map.begin();
    while( mapIt != m_Map.end() )
      {
      if( mapIt->second == 0 )
        { 
        // this value must be removed from the histogram
        // The value must be stored and the iterator updated before removing the value
        // or the iterator is invalidated.
        TInputPixel toErase = mapIt->first;
        mapIt++;
        m_Map.erase(toErase);
        }
      else
        {
        mapIt++;
        }
      }

    // and return the value
    if( !m_Map.empty() )
      {
      const TInputPixel & max = m_Map.rbegin()->first;
      const TInputPixel & min = m_Map.begin()->first;
      if( m_FlatToForeground )
        {
        if( centerPixel <= min )
          {
          return m_ForegroundValue;
          }
        }
      else
        {
        if( centerPixel <= min && centerPixel != max )
          {
          return m_ForegroundValue;
          }
        }
      }
    return m_BackgroundValue;
    }

  MapType m_Map;







  //
  // the vector based algorithm
  //

  inline void initVector()
    {
    // initialize members need for the vector based algorithm
    m_Vector.resize( static_cast<int>( NumericTraits< TInputPixel >::max() - NumericTraits< TInputPixel >::NonpositiveMin() + 1 ), 0 );
    m_Max = NumericTraits< TInputPixel >::NonpositiveMin();
    m_Min = NumericTraits< TInputPixel >::max();
    m_Count = 0;
    }

  inline void AddPixelVector( const TInputPixel &p )
    {
    m_Vector[ static_cast<int>( p - NumericTraits< TInputPixel >::NonpositiveMin() ) ]++;
    if( p > m_Max )
      { m_Max = p; }
    if( p < m_Min )
      { m_Min = p; }
    m_Count++;
    }

  inline void RemovePixelVector( const TInputPixel &p )
    {
    m_Vector[ static_cast<int>( p - NumericTraits< TInputPixel >::NonpositiveMin() ) ]--;
    m_Count--;
    if( m_Count > 0 )
      {
      while( m_Vector[ static_cast<int>( m_Max - NumericTraits< TInputPixel >::NonpositiveMin() ) ] == 0 )
        { m_Max--; }
      while( m_Vector[ static_cast<int>( m_Min - NumericTraits< TInputPixel >::NonpositiveMin() ) ] == 0 )
        { m_Min++; }
      }
    else
      {
      m_Max = NumericTraits< TInputPixel >::NonpositiveMin();
      m_Min = NumericTraits< TInputPixel >::max();
      }
    }

  inline const TOutputPixel & GetValueVector( const TInputPixel & centerPixel )
    {
    if( centerPixel > m_Threshold )
      {
      return m_BackgroundValue;
      }

    if( m_FlatToForeground )
      {
      if( m_Count > 0 &&  centerPixel <= m_Min )
        {
        return m_ForegroundValue;
        }
      }
    else
      {
      if( m_Count > 0 &&  centerPixel <= m_Min && centerPixel != m_Max )
        {
        return m_ForegroundValue;
        }
      }
      return m_BackgroundValue;
    }

  std::vector<unsigned long> m_Vector;
  TInputPixel m_Min;
  TInputPixel m_Max;
  unsigned long m_Count;


  void SetForegroundValue( TOutputPixel v )
    {
    m_ForegroundValue = v;
    }

  void SetBackgroundValue( TOutputPixel v )
    {
    m_BackgroundValue = v;
    }

  void SetFlatToForeground( bool v )
    {
    m_FlatToForeground = v;
    }

  void SetThreshold( TInputPixel v )
    {
    m_Threshold = v;
    }


  TOutputPixel m_ForegroundValue;
  TOutputPixel m_BackgroundValue;
  bool m_FlatToForeground;
  TInputPixel m_Threshold;

};
} // end namespace Function

/**
 * \class LocalMinimaImageFilter
 * \brief Find the pixels which are local extrema, according to a nighborhood
 *
 * 
 * \author Gaëtan Lehmann. Biologie du Développement et de la Reproduction, INRA de Jouy-en-Josas, France.
 *
 * \ingroup ImageEnhancement  MathematicalMorphologyImageFilters
 */


template<class TInputImage, class TOutputImage, class TKernel=Neighborhood< bool, TInputImage::ImageDimension > >
class ITK_EXPORT LocalMinimaImageFilter : 
    public MovingHistogramImageFilter<TInputImage, TOutputImage, TKernel,
      typename  Function::LocalMinimaHistogram< typename TInputImage::PixelType, typename TOutputImage::PixelType > >
{
public:
  /** Standard class typedefs. */
  typedef LocalMinimaImageFilter Self;
  typedef MovingHistogramImageFilter<TInputImage, TOutputImage, TKernel,
      typename  Function::LocalMinimaHistogram< typename TInputImage::PixelType, typename TOutputImage::PixelType > >  Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  /** Standard New method. */
  itkNewMacro(Self);  

  /** Runtime information support. */
  itkTypeMacro(LocalMinimaImageFilter, 
               ImageToImageFilter);
  
  /** Image related typedefs. */
  typedef TInputImage InputImageType;
  typedef TOutputImage OutputImageType;
  typedef typename TInputImage::RegionType RegionType ;
  typedef typename TInputImage::SizeType SizeType ;
  typedef typename TInputImage::IndexType IndexType ;
  typedef typename TInputImage::PixelType PixelType ;
  typedef typename TInputImage::OffsetType OffsetType ;
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
  typedef typename TOutputImage::PixelType OutputPixelType ;

  typedef typename Function::LocalMinimaHistogram< PixelType, OutputPixelType > HistogramType;


  /** Image related typedefs. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);


  /** Return true if the vector based algorithm is used, and
   * false if the map based algorithm is used */
  static bool GetUseVectorBasedAlgorithm()
    { return HistogramType::useVectorBasedAlgorithm(); }
  
  /**
   * Set/Get the value in the output image to consider as "foreground".
   * Defaults to maximum value of PixelType.
   */
  itkSetMacro(ForegroundValue, OutputPixelType);
  itkGetConstMacro(ForegroundValue, OutputPixelType);

  /**
   * Set/Get the value used as "background" in the output image.
   * Defaults to NumericTraits<PixelType>::NonpositiveMin().
   */
  itkSetMacro(BackgroundValue, OutputPixelType);
  itkGetConstMacro(BackgroundValue, OutputPixelType);

  itkSetMacro(FlatToForeground, bool);
  itkGetConstMacro(FlatToForeground, bool);
  itkBooleanMacro(FlatToForeground);

  itkSetMacro(Threshold, PixelType);
  itkGetConstMacro(Threshold, PixelType);

protected:
  LocalMinimaImageFilter()
    {
    m_ForegroundValue = NumericTraits<OutputPixelType>::max();
    m_BackgroundValue = NumericTraits<OutputPixelType>::NonpositiveMin();
    m_FlatToForeground = true;
    m_Threshold = NumericTraits<PixelType>::max();
    };

  ~LocalMinimaImageFilter() {};

  virtual HistogramType * NewHistogram()
    {
    HistogramType * histogram = Superclass::NewHistogram();
    histogram->SetForegroundValue( this->GetForegroundValue() );
    histogram->SetBackgroundValue( this->GetBackgroundValue() );
    histogram->SetFlatToForeground( this->GetFlatToForeground() );
    histogram->SetThreshold( this->GetThreshold() );
    return histogram;
    }

  void PrintSelf(std::ostream &os, Indent indent) const
    {
    Superclass::PrintSelf(os, indent);
    os << indent << "ForegroundValue: "  << static_cast<typename NumericTraits<OutputPixelType>::PrintType>(m_ForegroundValue) << std::endl;
    os << indent << "BackgroundValue: "  << static_cast<typename NumericTraits<OutputPixelType>::PrintType>(m_BackgroundValue) << std::endl;
    os << indent << "FlatToForeground: "  << m_FlatToForeground << std::endl;
    os << indent << "Threshold: "  << static_cast<typename NumericTraits<PixelType>::PrintType>(m_Threshold) << std::endl;
    }

private:
  LocalMinimaImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  OutputPixelType m_ForegroundValue;
  OutputPixelType m_BackgroundValue;
  bool m_FlatToForeground;
  PixelType m_Threshold;

} ; // end of class

} // end namespace itk
  
#endif


