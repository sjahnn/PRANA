## This code is to compute total connectivity of each gene from the association matrix.

## An input (asso.matinput) is association matrix.
## An output is total connectivity of each gene, denoted as \hat{\theta}_{k} in the paper.
thetahats = function(asso.matinput) {
  if(class(asso.matinput)[1] == "network") {
    asso.matinput <- get_association_matrix(asso.matinput)
  }
  results = vector()
  for(j in 1:ncol(asso.matinput)) {
    results[j] = sum(asso.matinput[j, -j])
  }
  return(results)
}


