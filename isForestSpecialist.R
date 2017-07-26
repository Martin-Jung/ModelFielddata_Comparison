isForestSpecialist <- function (inputID) 
{
  print(paste("Now processing ID",inputID))
  # Input data need to be joined first with IUCN ids
  # Queries the IUCN database if a species is associated with forests and has major importance = yes
  library(XML)
  library(RCurl)
  library(stringr)
  if(is.na(inputID)) return(NA) else 
    { 
      # Generate the link
      link <- paste("http://www.iucnredlist.org/details/classify/",inputID, "/0", sep = "")      
      # Download the html file
      webpage <- getURL(link)
      webpage <- readLines(tc <- textConnection(webpage)); close(tc)      
      # Parse site structure and extract body
      pagetree <- htmlTreeParse(webpage, error=function(...){})
      body <- pagetree$children$html$children$body
      # Then extract habitat container
      b <- try( xpathApply(body, "//div[@id='habitat']")[[1]] )
      if(class(b)[1] =="try-error") {
        print(paste("Species ID",inputID,"does not seem to have habitat info..."))
        return(NA)
      }
        
      # get text
      bb <- toString.XMLNode(b)
      # split at hr
      h2 <- unlist(str_split(bb,"<hr",n = Inf))
      # check if forest and mayor importance occur
      res <- vector()
      for(i in 1:length(h2)){        
        # Reduce does the trick 
        ix <- c( Reduce('&', lapply(c("Forest","Yes"), grepl, h2[i])),
           Reduce('&', lapply(c("forest","yes"), grepl, h2[i])),
           Reduce('&', lapply(c("Forest","yes"), grepl, h2[i])),
           Reduce('&', lapply(c("forest","Yes"), grepl, h2[i]))
        )
        if(any(ix)) res <- c(res,"Forest-specialist")
        iy <- c( Reduce('&', lapply(c("Forest","No"), grepl, h2[i])),
                 Reduce('&', lapply(c("forest","no"), grepl, h2[i])),
                 Reduce('&', lapply(c("Forest","no"), grepl, h2[i])),
                 Reduce('&', lapply(c("forest","No"), grepl, h2[i]))
        )
        if(any(iy)) res <- c(res,"Forest-associated")
    }
    # Check if one clear association. 
    if( "Forest-specialist" %in% res) return("Forest-specialist") else 
      if( length(unique(res)) == 2 ) return("Forest-specialist") else
        if( length(unique(res) == 1)) return("Forest-associated") else
          return("Other Habitats") # If none return NA by default if nothing else works    
  }
}