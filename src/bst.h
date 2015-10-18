/* 
 * Copyright (C) 2010    Niall Moran, Graham Kells, Jiri Vala
 * 
 * 
*/

#include "doqo.h"

#ifdef _BOOST_BST_
#include <boost/intrusive/avl_set.hpp>
#include <vector>
#include <algorithm>
#include <cassert>


//This is a base hook optimized for size
class CuPetscInt : public boost::intrusive::avl_set_base_hook<boost::intrusive::optimize_size<true> >
{
   public:
   uPetscInt int_;
   //This is a member hook
   boost::intrusive::avl_set_member_hook<> member_hook_;

   CuPetscInt() {}	
   CuPetscInt(uPetscInt i):int_(i){}
   void set_value(uPetscInt i ) {
   	int_ = i ;
   }
   friend bool operator< (const CuPetscInt &a, const CuPetscInt &b)
      {  return a.int_ < b.int_;  }
   friend bool operator> (const CuPetscInt &a, const CuPetscInt &b)
      {  return a.int_ > b.int_;  }
   friend bool operator== (const CuPetscInt &a, const CuPetscInt &b)
      {  return a.int_ == b.int_;  }
};

typedef boost::intrusive::avl_set<CuPetscInt>  BST;


#endif
