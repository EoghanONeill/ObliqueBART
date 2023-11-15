#include <Rcpp.h>
using namespace Rcpp;


// #include <Rcpp.h>
// #include <R_ext/BLAS.h>
// using namespace Rcpp;
//
// inline void mat_vec_mult_vanilla
//   (double const * __restrict__ m,
//    double const * __restrict__ v,
//    double * __restrict__ const res,
//    size_t const dn, size_t const dm) noexcept {
//   for(size_t j = 0; j < dm; ++j, ++v){
//     double * r = res;
//     for(size_t i = 0; i < dn; ++i, ++r, ++m)
//       *r += *m * *v;
//   }
// }
//
// inline void mat_vec_mult
//   (double const * __restrict__ const m,
//    double const * __restrict__ const v,
//    double * __restrict__ const res,
//    size_t const dn, size_t const dm) noexcept {
//   size_t j(0L);
//   double const * vj = v,
//     * mi = m;
//   constexpr size_t const ncl(8L);
//   {
//     double const * mvals[ncl];
//     size_t const end_j = dm - (dm % ncl),
//       inc = ncl * dn;
//     for(; j < end_j; j += ncl, vj += ncl, mi += inc){
//       double *r = res;
//       mvals[0] = mi;
//       for(size_t i = 1; i < ncl; ++i)
//         mvals[i] = mvals[i - 1L] + dn;
//       for(size_t i = 0; i < dn; ++i, ++r)
//         for(size_t ii = 0; ii < ncl; ++ii)
//           *r += *(vj + ii) * *mvals[ii]++;
//     }
//   }
//
//   mat_vec_mult_vanilla(mi, vj, res, dn, dm - j);
// }
//
// // [[Rcpp::export("mat_vec_mult", rng = false)]]
// NumericVector mat_vec_mult_cpp(NumericMatrix m, NumericVector v){
//   size_t const dn = m.nrow(),
//     dm = m.ncol();
//   NumericVector res(dn);
//   mat_vec_mult(&m[0], &v[0], &res[0], dn, dm);
//   return res;
// }
//
// // [[Rcpp::export("mat_vec_mult_vanilla", rng = false)]]
// NumericVector mat_vec_mult_vanilla_cpp(NumericMatrix m, NumericVector v){
//   size_t const dn = m.nrow(),
//     dm = m.ncol();
//   NumericVector res(dn);
//   mat_vec_mult_vanilla(&m[0], &v[0], &res[0], dn, dm);
//   return res;
// }
//
// // [[Rcpp::export(rng = false)]]
// NumericVector my_mm(NumericMatrix m, NumericVector v){
//   int nRow = m.rows();
//   int nCol = m.cols();
//   NumericVector ans(nRow);
//   double v_j;
//   for(int j = 0; j < nCol; j++){
//     v_j = v[j];
//     for(int i = 0; i < nRow; i++){
//       ans[i] += m(i,j) * v_j;
//     }
//   }
//   return(ans);
// }
//
// // [[Rcpp::export(rng = false)]]
// NumericVector blas_mm(NumericMatrix m, NumericVector v){
//   int nRow = m.rows();
//   int nCol = m.cols();
//   NumericVector ans(nRow);
//   char trans = 'N';
//   double one = 1.0, zero = 0.0;
//   int ione = 1;
//   F77_CALL(dgemv)(&trans, &nRow, &nCol, &one, m.begin(), &nRow, v.begin(),
//            &ione, &zero, ans.begin(), &ione);
//   return ans;
// }
