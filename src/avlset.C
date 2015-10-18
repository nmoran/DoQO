#include <boost/intrusive/avl_set.hpp>
#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <list>
#include <ctime>
#include <cstdlib>
#include <ctime>


using namespace boost::intrusive;

//This is a base hook optimized for size
class MyClass : public avl_set_base_hook<optimize_size<true> >
{
   public:
   unsigned int int_;
   
   //This is a member hook
   avl_set_member_hook<> member_hook_;

   MyClass(){}

   MyClass(unsigned int i)
      :  int_(i)
      {}
      
   void set_value(unsigned int i) {
   	int_ = i;
   } 
   friend bool operator< (const MyClass &a, const MyClass &b)
      {  return a.int_ < b.int_;  }
   friend bool operator> (const MyClass &a, const MyClass &b)
      {  return a.int_ > b.int_;  }
   friend bool operator== (const MyClass &a, const MyClass &b)
      {  return a.int_ == b.int_;  }
};

//Define an avl_set using the base hook that will store values in reverse order
//typedef avl_set< MyClass, compare<std::greater<MyClass> > >     BaseSet;
typedef avl_set< MyClass,constant_time_size<true>, size_type<float> >  BaseSet;

using namespace std;

int insert_into_sorted_list(std::list<int> *sorted_list, int  item) {
	//PetscPrintf(PETSC_COMM_WORLD,"empty %d, size %d, true %d\n", sorted_list->empty(), sorted_list->size(),true);
	if ( sorted_list->empty() ) { 
		//PetscPrintf(PETSC_COMM_WORLD,"List address %d, item %d\n",sorted_list,item);
		sorted_list->push_back(item);
		//PetscPrintf(PETSC_COMM_WORLD,"List size %d, empty %d, first %d\n", sorted_list->size(), sorted_list->empty(), sorted_list->front());
		return 1;
	}
	std::list<int>::iterator iter = sorted_list->begin();
	while ( iter != sorted_list->end() && item > *iter ) {
		iter++;
	}
	if (iter != sorted_list->end() && item == *iter) {
		return 0;
	} else {
		//PetscPrintf(PETSC_COMM_WORLD,"Inserting %d\n",item);
		sorted_list->insert(iter,item);
		return 1;
	}
}


int main()
{
	int size = 1000;
	srand((unsigned)time(0)); 
   	BaseSet **baseset;
   	
   	//baseset = new BaseSet*[10];
   	baseset = (BaseSet**) malloc(sizeof(BaseSet*)*10);
   	   	
   	//baseset[5] = new BaseSet;
  	baseset[5] = (BaseSet*)malloc(sizeof(BaseSet));
  
  	std::vector<BaseSet> vec;
  	
  	//BaseSet bar*;
  	
  	//bar = new Baseset;
  	
  	//vec.push_back(bar);
  
	std::list<int > int_list ;
	
	clock_t start,finish;
	start = clock();
	for ( int  i = 0 ; i < size ; i++ ) {
	  insert_into_sorted_list(&int_list,rand()%(size*10));
	}  
	finish = clock();
	cout << "With list took : " <<  (double(finish)-double(start))/CLOCKS_PER_SEC << endl;
	
	MyClass *foo;
	int count = 0;
	start = clock();
	
	foo = (MyClass*)malloc(sizeof(MyClass));
	foo->set_value(0);
	//cout << foo->int_ << std::endl; 
	if  ( baseset[5]->insert(*foo).second == false ) {
		//delete foo;
		free(foo);
	} else {
		count ++ ;
	}
	
	for ( int i = 0 ; i < size ; i++ ) {
		//foo = new MyClass(rand()%(size*10));
		foo = (MyClass*)malloc(sizeof(MyClass));
		foo->set_value(rand()%(size*10));
		//cout << foo->int_ << std::endl; 
		if  ( baseset[5]->insert(*foo).second == false ) {
			//delete foo;
			free(foo);
		} else {
			count ++ ;
		}
	}  
	finish = clock();
	std::cout << "With tree took : "<<  (double(finish)-double(start))/CLOCKS_PER_SEC << std::endl;
	
	cout << "Tree size: " << baseset[5]->size() << " , " << count << std::endl;
	
	BaseSet::iterator it(baseset[5]->begin());
	//for ( ; it != baseset.end() ; it++ ) {
	while ( it != baseset[5]->end() ) {
		//std::cout << (*it).int_ << std::endl;
		MyClass *element = &(*it);
		it = baseset[5]->erase(it);
		free(element);
	}
	
	return 0;
}