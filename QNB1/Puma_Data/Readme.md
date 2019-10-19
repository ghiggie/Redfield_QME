The data contained in this folder was generated with a fixed
Simpson's Rule. In bath.f95, I had an issue in which the counter
vairable M was not always even. To fix this, I check to see
whether M is even before proceeding. If it is, then proceed to
the loop. If it is not, increment M by 1 and re-evaluate dt using
dt = t / M. This dt will be different than the dt given in the 
function definition, but it will be close enough that the accuracy
will not be affected to any appreciable degree. Once dt and M are
redefined, proceed with the loop.

In Run1, the coefficients of the bath correlation function are
unaltered.

In Run2, the coefficients c1 and c2 (corresponding to the
exponential terms) are halved.

In Run3, all the coefficients are halved.
