#include <list>
using namespace std;

template<class ______>
typename ______::value_type _(______ &____)
{
  if (____.size()==0)
    return 0;
  else
    {
      typename ______:: iterator ___=____.begin();
      typename ______::value_type __=____.front();
      ____.erase(___);
      return __+_(____);
    }
    }
