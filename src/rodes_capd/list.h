/*   File: mult_list.h 
 
     A template version of a linked list, 
     with dynamic memory allocation.

     Usage: List<int>, List<vector>, etc...
 
     Latest edit: Sun Jul 24 2000
*/

#ifndef MULT_LIST_H
#define MULT_LIST_H 

#include <iostream.h>
#include <stdio.h>

//#include "list_funcs.h"

////////////////////////////////////////////////////////////////////////////
// Prototype class and friend functions
template < class TYPE > class List;
template < class TYPE > TYPE & First ( List < TYPE > & );
template < class TYPE > TYPE & Next ( List < TYPE > & );
template < class TYPE > TYPE & Last ( const List < TYPE > & ); 
template < class TYPE > TYPE & Current ( const List < TYPE > & );
template < class TYPE > void RemoveCurrent ( List < TYPE > & );
template < class TYPE > 
ostream & operator << ( ostream &, const List < TYPE > & );

////////////////////////////////////////////////////////////////////////////

template<class TYPE>
struct ListObject 
{
  ListObject  *next;
  TYPE  element;
};

////////////////////////////////////////////////////////////////////////////
// 'friend' members are there so that they have access to speific
// private members of List
template<class TYPE>
class List 
{
  //private:
  public:
    ListObject<TYPE> *start;	 
    ListObject<TYPE> *end;	 
    ListObject<TYPE> *current;	 
    ListObject<TYPE> *lastcur;	 
    int  len;		 
    int  maxlen;		 
    //public:
  // initialize an empty list
    List()
      { start = end = current = lastcur = __null ; len = maxlen = 0; }
    ~List();
    List(List & li)
      { printf("List (Copy): only references allowed\n"); exit(1); }
    friend int Finished     (const List & li)
      { return (li.current == __null); }
    friend int IsEmpty      (const List & li)
      { return (li.start == __null ); }
    friend int Length       (const List & li)
      { return li.len; }
    friend int MaxLength    (const List & li)
      { return li.maxlen; }
    friend void ResetLength (List & li)
      { li.maxlen = li.len; }
    void operator +=  (const TYPE &);
    void operator *=  (const TYPE &);
    void operator --  ();

    friend TYPE & First   <TYPE>(List &);
    friend TYPE & Next    <TYPE>(List &);
    friend TYPE & Last    <TYPE>(const List &);
    friend TYPE & Current <TYPE>(const List &);
    friend void RemoveCurrent <TYPE>(List &);
    friend ostream & operator << <TYPE>(ostream &, const List &);

};

////////////////////////////////////////////////////////////////////////////
// Member function definitions follow

template<class TYPE>
List<TYPE>::~List()		 
{
  ListObject<TYPE> *temp;

  while (start != __null) {
    temp = start; start = start->next;
    delete temp;
  }
  start = end = __null;
  len = 0;
}

////////////////////////////////////////////////////////////////////////////

template<class TYPE>
void List<TYPE>::operator += (const TYPE & obj)
{
  ListObject<TYPE> *p = new ListObject<TYPE>;
  p->element = obj;
  p->next = __null;
  if (start == __null) start = p;
  else end->next = p;
  end = p;
  len++;
  if (len > maxlen) maxlen = len;
}

////////////////////////////////////////////////////////////////////////////

template<class TYPE>
void List<TYPE>::operator *= (const TYPE & obj)
{
  ListObject<TYPE> *p = new ListObject<TYPE>;
  p->element = obj;
  p->next = start;
  if (start == __null) end = p;
  start = p;
  len++;
  if (len > maxlen) maxlen = len;
}

////////////////////////////////////////////////////////////////////////////

template<class TYPE>
void List<TYPE>::operator -- ()
{
  if (start == __null) printf("List (--): empty list\n");
  ListObject<TYPE> *p = start->next;
  delete start;
  start = p;
  if (start == __null) end = __null;
  len--;
}

////////////////////////////////////////////////////////////////////////////

template<class TYPE>
ostream & operator << (ostream & o, const List<TYPE> &li)
{
  int  i = 0;
  if (li.start == __null) o << "*EMPTY*\n";
  else {
    for (ListObject<TYPE> *p = li.start; p != __null; p = p->next) {
      o.width(3);
      o << (++i) << ": ";
      o << (p->element) << '\n';
    }
  }
  return o;
}

////////////////////////////////////////////////////////////////////////////

template<class TYPE>
TYPE & First (List<TYPE> & li)
{
  li.current = li.start;
  li.lastcur = __null;
  if (li.current == __null) printf("List (First): empty list\n");
  return li.current->element;
}

////////////////////////////////////////////////////////////////////////////

template<class TYPE>
TYPE & Next (List<TYPE> & li)
{
  if ((li.current == __null) && (li.lastcur == __null)) return First(li);
  if (li.current == __null) return li.lastcur->element;
  li.lastcur = li.current;
  li.current = li.current->next;
  return (li.current == __null) ? li.lastcur->element : li.current->element;
}

////////////////////////////////////////////////////////////////////////////

template<class TYPE>
TYPE & Last (const List<TYPE> & li) 
{
  if (li.end == __null) printf("List (Last): empty list\n");
  return li.end->element;
}

////////////////////////////////////////////////////////////////////////////

template<class TYPE>
TYPE & Current (const List<TYPE> & li) 
{
  if (li.current == __null) printf("List (Current): no element\n");
  return li.current->element;
}

////////////////////////////////////////////////////////////////////////////

template<class TYPE>
void RemoveCurrent (List<TYPE> & li)
{
  ListObject<TYPE> *del_cur = li.current;
  if (li.current == __null) printf("List (RemoveCurrent): no element\n");
  li.current = li.current->next;
  delete del_cur;
  li.len--;
  if (li.lastcur == __null) li.start = li.current;
  else li.lastcur->next = li.current;
  if (li.current == __null) li.end = li.lastcur;
}

//#endif // __linux

////////////////////////////////////////////////////////////////////////////

#endif // MULT_LIST_H 
