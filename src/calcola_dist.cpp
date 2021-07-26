#include <Rcpp.h>
#include <math.h>       /* round, floor, ceil, trunc */

using namespace Rcpp;
// Below is a simple example of exporting a C++ function to R.
// You can source this function into an R session using the
// Rcpp::sourceCpp function (or via the Source button on the
// editor toolbar)
// For more on using Rcpp click the Help button on the editor
// toolbar



// [[Rcpp::export]]
NumericVector cal_dist(NumericVector x1,NumericVector x2,int ncont, int nbin){

	NumericVector out(3,0.0);
	//Rcpp::Rcout<<"out="<<out<<"\n";

	int j;
	for(j=0;j<ncont;j++){
		out[1] = out[1]+  pow(x1[j]-x2[j],2.);
	}
	//std::cout<<"ciao!!"<<__LINE__<<std::endl;
	out[1] =std::pow(out[1],0.5);

	for(j=(ncont);j<(ncont+nbin);j++){
		out[2] = out[2]+ (x1[j]!=x2[j]);
	}

	out[0]=out[1]+out[2];

	//std::cout<<"ciao!!"<<__LINE__<<std::endl;
	
	return(out);

}






// [[Rcpp::export]]
List dist_r(NumericMatrix xx,int ncont, int nbin) {

	int n= xx.nrow();
	int p= xx.ncol();


	if(p!=ncont+nbin){
		Rcpp::stop("errore numero di binarie +numero di continue non uguale a numero di regressori");
	}
	int M = (std::pow(n,2)-n)/2;
	std::cout<<"CIAO! Il numero M"<<n<<" "<<p<<" "<<M<<std::endl;
	NumericVector distanze(M);
	NumericVector euclidea(M);
	NumericVector hamming(M);

	NumericVector appoggio(3);
	int iterator=0;
	NumericVector xi(p);
	NumericVector xj(p);
	//std::cout<<"RI CIAO!"<<n<<" "<<p<<" "<<M<<std::endl;

	Rcpp::IntegerVector I(M);
	Rcpp::IntegerVector J(M);
	for(int i =0;i<n;i++){
		for(int j=(i+1);j<n;j++){
			//std::cout<<"ciao!!"<<iterator<<"; i=" <<i<<" j="<<j  <<std::endl;
			xi =xx(i,_);
			xj= xx(j,_);
			appoggio=cal_dist(xi,xj,ncont,nbin);
			//Rcpp::Rcout<<"appoggio="<<appoggio <<std::endl;

			distanze[iterator]=appoggio[0];
			euclidea[iterator]=appoggio[1];
			hamming[iterator]=appoggio[2];
			I[iterator]=i+1;
			J[iterator]=j+1;
			iterator +=1;
		}
		std::cout<<"CIAO i="<<i<<std::endl;

	}

// 	Rcpp::List ret;
// 	ret["distanze"] = distanze;
// 	ret["I"] =I;
// 	//ret["J"]=J;
// 	return ret;
// 

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
		
                          //Rcpp::Named("lst", I),
                          //Rcpp::Named("vec2", xj));


}



// [[Rcpp::export]]

NumericVector centroide(NumericMatrix xx,int ncont, int nbin){
	int n= xx.nrow();
	int p= xx.ncol();
	if(p!=ncont+nbin){
		Rcpp::stop("errore numero di binarie +numero di continue non uguale a numero di regressori");
	}
	Rcpp::NumericVector out(p);
        //Rcpp::Rcout << "The value out is " << out << std::endl;
	for(int i=0;i<n;i++){
	out += xx(i,_);
	}
	
	out=out/n;
	for(int j=(ncont);j<(ncont+nbin);j++){
	out[j]=round(out[j]);
	}
	return(out);
}

// [[Rcpp::export]]
double calcola_D(NumericMatrix xx,int ncont, int nbin){
	NumericVector centro=centroide(xx,ncont,nbin);
	double out=0;
	int n= xx.nrow();

	NumericVector pippo;
	for(int i=0;i<n;i++){
		pippo= cal_dist(centro,xx(i,_),ncont,nbin);
		out += pippo[0];
	}
	return(out);



}

// [[Rcpp::export]]
double calcola_D_norm(NumericMatrix xx,int ncont, int nbin){
	NumericVector centro=centroide(xx,ncont,nbin);
	double out=0;
	int n= xx.nrow();

	NumericVector pippo;
	for(int i=0;i<n;i++){
		pippo= cal_dist(centro,xx(i,_),ncont,nbin);
		out+= pippo[0];
	}
	return(pow(n,-1)*out);



}

// 
// 
// 
// 
// // Copy the second column into new object
// //NumericVector zz1 = xx( _, 1);
// // Fill with value
// //int xsize = xx.nrow() * xx.ncol();
// 
// 
// 
// 
// 
// 
// //  return Rcpp::List::create(Rcpp::Named("A") = A,
// //                             Rcpp::Named("PH") = PH,
// //                             Rcpp::Named("Z") = Z,
// //                             Rcpp::Named("Zi") = Zi,
// //                             Rcpp::Named("Ze") = Ze);
// // 
