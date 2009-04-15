## Normalize distance using the following codes (same as in C code)
## 0 => no normalization
## 1 => Abbott
## 2 => Elzinga
## 3 => Maximum possible distance normalization (divide by maxdist)



normdist <- function(rawdist, maxdist, l1, l2, norm) {

	if (rawdist==0) {
		return(0)
 	}
	if (norm==0) { #Without normalization
			return(rawdist)
	} else if (norm==1) { #Abbott normalization
			if (l1>l2) return(rawdist/l1)
			else if (l2>0) return(rawdist/l2)
			return(0)
	} else if (norm==2) { #Elzinga normalization
			if (l1*l2==0) {
				if (l1!=l2) return(1)
				return(0)
			}
			return(1-((maxdist-rawdist)/(2*sqrt(l1*l2))))
	} else if (norm==3) { #Maximum possible distance normalization
			if (maxdist==0) return(1)
			return(rawdist/maxdist)
	}
	stop("Unknow distance normalization used")
}