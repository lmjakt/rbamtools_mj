% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% as.data.frame
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\name{as.data.frame-methods}
\docType{methods}
\alias{as.data.frame-methods}
\title{Conversion of bamRange or gapList into a data.frame }
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{
\code{as.data.frame} functions convert objects of type bamRange or gapList
into a data.frame.
Data from each element (bamAlign or Align-gap) is written into a single line.
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Methods
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\section{Methods}{
\describe{
    \item{\code{signature(x = "bamRange")}}{
        Conversion of bamRange into data.frame}
    %
    \item{\code{signature(x = "gapList" )}}{
        Conversion of gapList into data.frame}
}}
\keyword{methods}
\keyword{bamRange}
\keyword{bamAlign}
\keyword{gapList}
