// Copyright 2008-2016 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup mtSpGlue
//! @{



template<typename out_eT, typename T1, typename T2, typename spglue_type>
class mtSpGlue : public SpBase<out_eT, mtSpGlue<out_eT, T1, T2, spglue_type> >
  {
  public:
  
  typedef          out_eT                       elem_type;
  typedef typename get_pod_type<out_eT>::result pod_type;
  
  static const bool is_row  = spglue_type::template traits<T1,T2>::is_row;
  static const bool is_col  = spglue_type::template traits<T1,T2>::is_col;
  static const bool is_xvec = spglue_type::template traits<T1,T2>::is_xvec;
  
  inline  mtSpGlue(const T1& in_A, const T2& in_B);
  inline ~mtSpGlue();
  
  template<typename eT2>
  arma_inline bool is_alias(const SpMat<eT2>& X) const;
  
  const T1& A;    //!< first operand;  must be derived from SpBase
  const T2& B;    //!< second operand; must be derived from SpBase
  };



//! @}