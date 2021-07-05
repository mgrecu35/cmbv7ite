#include "TKheaders.h"

L2AKu_NS        kuswath;
L2AKuX_FS        kuxswath;

void read1scan()
{
  extern TKINFO granuleHandle2AKu;
  extern TKINFO tkfileinfo;
  extern TKINFO dprtkfile;
  int status;
  extern L2ADPR_SWATHS dprswath;
  extern L2AKu_NS      kuswath;
  while (TKendOfFile (&granuleHandle2AKu) != TK_EOF)
    {
      status=TKreadScan(&granuleHandle2AKu,&kuswath);
      status = TKreadScan(&dprtkfile,&dprswath);
    }
}

void read1scanx()
{
  extern TKINFO granuleHandle2AKu;
  extern TKINFO tkfileinfo;
  extern TKINFO dprtkfile;
  int status;
  extern L2ADPRX_SWATHS dprxswath;
  extern L2AKuX_FS      kuxswath;
  while (TKendOfFile (&granuleHandle2AKu) != TK_EOF)
    {
      status=TKreadScan(&granuleHandle2AKu,&kuxswath);
      status = TKreadScan(&dprtkfile,&dprxswath);
    }
}
