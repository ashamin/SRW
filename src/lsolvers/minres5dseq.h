#ifndef MINRES_5D_SEQ_H
#define MINRES_5D_SEQ_H

#include "minres5d.h"

class minres5dseq : public minres5d{
public:
  minres5dseq();
  double* solve();
};

#endif