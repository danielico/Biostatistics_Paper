
# Modèles simulations 
#--------------------------------------------------------------

list.model.simu <- function( knot.agec = 0 )
{
  ms.lin.ph   = "~ unsurt + logt + intnum + I(intnum^2) + agec"
  ms.lin.nph  = "~ unsurt + logt + intnum + I(intnum^2) + agec + agec:intnum + agec:I(intnum^2) + agec:I(intnum^3) + agec:I((intnum-1)^3*(intnum>1))"
  ms.nlin.ph  = paste( "~ unsurt + logt + intnum + I(intnum^2) +" , "agec + I(agec^2) + I(agec^3) + I((agec-",knot.agec,")^3*(agec>",knot.agec,"))")
  
  return( list( ms.lin.ph = ms.lin.ph , ms.lin.nph = ms.lin.nph , ms.nlin.ph = ms.nlin.ph ) )
}


# Modeles analyses
#------------------------------------------------

list.model.ana <- function( knot.agec = 0 )
{
  ma.lin.ph   = "~ intnum + I(intnum^2) + I((intnum-1)^2*(intum>1)) + I((intnum-5)^2*(intum>5)) + agec"
  ma.lin.nph  = "~ intnum + I(intnum^2) + I((intnum-1)^2*(intum>1)) + I((intnum-5)^2*(intum>5)) + agec:intnum + agec:I(intnum^2) + agec:I(intnum^3) + agec:I((intnum-1)^3*(intnum>1))"
  ma.nlin.ph  = paste( "~ intnum + I(intnum^2) + I((intnum-1)^2*(intum>1)) + I((intnum-5)^2*(intum>5))" , "+ agec + I(agec^2) + I(agec^3) + I((agec-",knot.agec,")^3*(agec>",knot.agec,"))")
   
  return( list( ma.lin.ph = ma.lin.ph , ma.lin.nph = ma.lin.nph , ma.nlin.ph = ma.nlin.ph ))
}
