\name{logfam}
\alias{logfam}
\title{Calculated the logarithm of the Marginal Likelihood of Gaussian multivariate regression under the DAG model}
\usage{
logfam(Y,D,X)
}
\description{
Calculated the logarithm of the Marginal Likelihood of Gaussian multivariate regression under the DAG model.
}
\arguments{
\item{Y}{ response matrix \eqn{n \times q} that regressed on the \eqn{X} matrix represents a Gaussian multivariate variable}
\item{D}{ adjacency matrix \eqn{q \times q} corresponds to the DAG structure. It should have \eqn{d_{ii}=0}; if the node \eqn{j} points to the node \eqn{i} then \eqn{d_{ij}=1} and \eqn{d_{ji}=0}; if there is not a link between nodes \eqn{i} and \eqn{j} then \eqn{d_{ij}=d_{ji}=0}}
\item{X}{ regression matrix \eqn{p \times q} with \eqn{p} predictors}
}

\author{
Cecilia Piergiovanna Trevisi c.trevisi1@campus.unimib.it
}

\details{The expression of the marginal density of the matrix \eqn{Y} under the DAG \eqn{D}:\cr{}
\eqn{
\textit{m}_\textit{D}(Y)=\prod_{j=1}^{q} \frac{m(Y_{faD(j)})}{m(Y_{paD(j)}}
}\cr{}
where \eqn{D} is the DAG formed by \eqn{
  \{1, \ldots, q\}
} vertices; those pointing to general vertex \eqn{j} are called parents (denoted by \eqn{
  pa_D(j)}), and \eqn{y_{{ipa}_D}(j)} is the subvector of \eqn{y_i} indexed by \eqn{pa_D(j)}.\cr{} The \eqn{fa_D(j)=pa_D(j) \cup \{j\}} is the family of \eqn{j} in \eqn{D}.\cr{}
\eqn{m(Y_{faD(j)})} and \eqn{m(Y_{paD(j)})} are calculated using the following formula:\cr{}
\eqn{
\textit{m}(Y_J)= \pi^{-\frac{(n-n_0)|J|}{2}} \frac{\Gamma_{|J|}\left(\frac{a_D+n-p-1-|\overline{J}|}{2}\right)}{\Gamma_{|J|}\left(\frac{a_D+n_0-p-1-|\overline{J}|}{2}\right)} \left(\frac{n_0}{n}\right)^{\frac{|J|(a_D+n_0-\overline{J})}{2}}|\hat{E_J}^T\hat{E_J}|^{-\frac{n-n_0}{2}}
}\cr{}
which represents the marginal distribution of the \eqn{Y_J} sub-matrix.\cr{}
\eqn{n_0} has been taken as \eqn{p+2} (p predictors) and \eqn{a_D} as \eqn{q-1} (q response variables).\cr{}
In order to simplify the calculations and better interpret the results, the expression of the marginal probability of the matrix Y under the DAG \eqn{D} has been applied to the logarithm.

}


\references{
Consonni, G. and La Rocca, L. and Peluso,S.(2016), Objective Bayes Covariate-Adjusted Sparse Graphical Model Selection
}

\examples{
install.packages('MASS')
install.packages('BDgraph')
library(BDgraph)
library(MASS)
n<-6  # 6 observations
p<-2  # 2 predictors
q<-3  # 3 response variables
X<-matrix(c(1,1,2,3,4,5,1,2,3,1,1,1),n,p)
X=cbind(X,1)  # regression matrix (2 predictors plus intercept)
B<-matrix(c(1,3,2,1,1,1,3,4,5),p+1,q) # coefficient matrix
M=\code{X \%*\% B} # mean matrix
Dw<-matrix(c(0,0,0,1,0,0,1,1,0),q,q) # matrix of complete graph
set.seed(123)
omega=rgwish( n = 1, adj = Dw, b = 10, D = diag(3)) # precision matrix
sigma=solve(omega)  # covariance matrix
set.seed(123)
Y<-t(apply(M,1,function(r) mvrnorm(1,r,sigma))) # response matrix

D<-matrix(c(0,1,0,0,0,0,1,1,0),q,q) # DAG matrix

X <- X[, 1:(ncol(X) - 1)] # remove last column of matrix X (consisting of only 1)
                        # because it is added during the calculation
                        # of the log-marginal likelihood

logfam(Y,D,X)
}
