/*
 * The code in this file is by Raffaele Argiento. I have only commented some 
 * functions.
 */

#include <Rcpp.h>
#include <math.h>       /* round, floor, ceil, trunc */

using namespace Rcpp;

//' Computes distances between vectors
//' 
//' @param x1 first vector
//' @param x2 second vector
//' @param ncont number of continuous covariates in vector x1 (or x2)
//' @param ncont number of binary covariates in vector x1 (or x2)
//' @return vector of length 3, it contains three distances that are 
//' sum of Euclidean and Hamming distance, Euclidean and Hamming distance, 
//' respectively.
// [[Rcpp::export]]
NumericVector cal_dist(NumericVector x1, NumericVector x2, int ncont, int nbin){

	NumericVector out(3,0.0);

	int j;
	for(j = 0; j < ncont; j++){
		out[1] = out[1] + pow(x1[j]-x2[j],2.);
	}
	out[1] =std::pow(out[1],0.5);

	for(j=(ncont);j<(ncont+nbin);j++){
		out[2] = out[2]+ (x1[j]!=x2[j]);
	}

	out[0]=out[1]+out[2];
	
	return(out);

}

//' Computes distance among points
//'
//' @param xx cluster covariates
//' @param ncont number of continuous covariates
//' @param nbin number of binary covariates
//' @return list of length 5. First three elements are sum of Euclidean and 
//' Hamming distance, Euclidean, Hamming distance respectively. Last two 
//' elements are vectors reporting indices of the elements whose distance is
//' computed.
// [[Rcpp::export]]
List dist_r(NumericMatrix xx, int ncont, int nbin) {

	int n = xx.nrow();
	int p = xx.ncol();


	if(p!=ncont+nbin){
		Rcpp::stop("errore numero di binarie +numero di continue non uguale a numero di regressori");
	}
	
	//number of distances to be computed
	int M = (std::pow(n,2)-n)/2;
	//std::cout<<"CIAO! Il numero M"<<n<<" "<<p<<" "<<M<<std::endl;
	NumericVector distanze(M);
	NumericVector euclidea(M);
	NumericVector hamming(M);

	NumericVector appoggio(3);
	int iterator=0;
	NumericVector xi(p);
	NumericVector xj(p);

	Rcpp::IntegerVector I(M);
	Rcpp::IntegerVector J(M);
	for(int i = 0; i < n; i++){
		for(int j = (i+1); j < n; j++){
			xi =xx(i,_);
			xj= xx(j,_);
			appoggio = cal_dist(xi,xj,ncont,nbin);

			distanze[iterator] = appoggio[0];
			euclidea[iterator] = appoggio[1];
			hamming[iterator] = appoggio[2];
			I[iterator] = i+1;
			J[iterator] = j+1;
			iterator += 1;
			}
		}
	
	List A=Rcpp::List::create(
    Rcpp::Named("dist")=NumericVector::create(NA_REAL),
    Rcpp::Named("euclidea")=NumericVector::create(NA_REAL),
    Rcpp::Named("hamming")=NumericVector::create(NA_REAL),
    Rcpp::Named("I")=IntegerVector::create(NA_INTEGER),
    Rcpp::Named("J") = IntegerVector::create(NA_INTEGER));
	
	A["dist"]=distanze;
	A["euclidea"]=euclidea;
	A["hamming"]=hamming;
	A["I"]=I;
	A["J"]=J;
	
	return A;
	}


//' Computes the centroid of a cluster
//'
//' @param xx cluster covariates
//' @param ncont number of continuous covariates
//' @param nbin number of binary covariates
//' @return a vector of lenght ncont+nbin that is the centroid. it is the 
//' marginal mean for continuous covariates and the median for binary ones
// [[Rcpp::export]]
NumericVector centroide(NumericMatrix xx, int ncont, int nbin){
	int n= xx.nrow();
	int p= xx.ncol();
	
	if(p!=ncont+nbin){
		Rcpp::stop("errore numero di binarie +numero di continue non uguale a numero di regressori");
	}
	Rcpp::NumericVector out(p);

	for(int i = 0; i < n; i++){
	  out += xx(i,_);
	}
	
	out = out/n;
	for(int j = (ncont); j < (ncont+nbin); j++){
	  out[j] = round(out[j]);
	}
	return(out);
}

//' Computes the sum of a cluster distance from its centroid
//'
//' @param xx cluster covariates
//' @param ncont number of continuous covariates
//' @param nbin number of binary covariates
//' @return a double that is the sum of the distances of each point in the 
//' cluster from the centroid of the cluster
// [[Rcpp::export]]
double calcola_D(NumericMatrix xx, int ncont, int nbin){
	NumericVector centro = centroide(xx,ncont,nbin);
	double out = 0;
	int n = xx.nrow();

	NumericVector pippo;
	for(int i = 0; i < n; i++){
		pippo = cal_dist(centro, xx(i,_), ncont, nbin);
		out += pippo[0];
	}
	return(out);
}

//' Computes the normalized sum of a cluster distance from its centroid
//'
//' @param xx cluster covariates
//' @param ncont number of continuous covariates
//' @param nbin number of binary covariates
//' @return a double that is the normalized sum of the distances of each point 
//' in the cluster from the centroid of the cluster
// [[Rcpp::export]]
double calcola_D_norm(NumericMatrix xx, int ncont, int nbin){
	NumericVector centro = centroide(xx,ncont,nbin);
	double out = 0;
	int n = xx.nrow();

	NumericVector pippo;
	for(int i = 0; i<n; i++){
		pippo = cal_dist(centro, xx(i,_), ncont, nbin);
		out += pippo[0];
	}
	return(pow(n,-1)*out);
}
