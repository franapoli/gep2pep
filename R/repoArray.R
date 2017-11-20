
### =========================================================================
### HDF5Array objects
### -------------------------------------------------------------------------


setClass("repoArraySeed",
         contains="array",
         representation(
             repopath="character", 
             name="character",
             dim="numeric"
         )
)

.validate_repoArraySeed <- function(x)
{
    ## ## 'rle' slot.
    ## if (!is(x@repo, "repo"))
    ##     return(wmsg("'x@repo' must be a repo object"))
    ## ## 'dim' slot.
    ## if (!is.integer(x@dim))
    ##     return(wmsg("'x@dim' must be an integer vector"))
    ## x_ndim <- length(x@dim)
    ## if (x_ndim == 0L)
    ##     return(wmsg("'x@dim' cannot be empty"))
    ## if (S4Vectors:::anyMissingOrOutside(x@dim, 0L))
    ##     return(wmsg("'x@dim' cannot contain negative or NA values"))

    TRUE
}

setValidity2("repoArraySeed", .validate_repoArraySeed)

setMethod("dim", "repoArraySeed", function(x) x@dim)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### subset_seed_as_array()
###

.subset_repoArraySeed_as_array <- function(seed, index)
{
##    return(integer())
    ##  print(index)
    print(1)
    rp <- repo_open(seed@repopath)
    print(2)
    x <- rp$get(seed@name)
    print(3)
    seed@dim <- 1:2
    print(4)
    ##dim(x) <- 2
    print(5)
    x <- rbind(x$ES,x$PV)
    print(dim(x))
    ans <- integer(4)
    dim(ans) <- 4
    return(ans)
}

setMethod("subset_seed_as_array", "repoArraySeed",
    .subset_repoArraySeed_as_array
)


### Return a HDF5ArraySeed object with NO dimnames!
### FIXME: Investigate the possiblity to store the dimnames in the HDF5 file
### and make dimnames() on the object returned by HDF5ArraySeed() bring them
### back.
repoArraySeed <- function(repopath, name)
{
    new2("repoArraySeed",
         repopath=repopath,
         name=name,
         dim=1)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### HDF5Array and HDF5Matrix objects
###
### We define these classes only for cosmetic reasons i.e. to hide the
### DelayedArray and DelayedMatrix classes from the user. The user will see
### and manipulate HDF5Array and HDF5Matrix objects instead of DelayedArray
### and DelayedMatrix objects.
###

setClass("repoArray", contains="DelayedArray")

setClass("repoMatrix", contains=c("DelayedMatrix", "repoArray"))

### Automatic coercion method from HDF5Array to HDF5Matrix silently returns
### a broken object (unfortunately these dummy automatic coercion methods don't
### bother to validate the object they return). So we overwrite it.
setAs("repoArray", "repoMatrix", function(from) new("repoMatrix", from))


setAs("ANY", "repoMatrix",
    function(from) as(as(from, "repoArray"), "repoMatrix")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

setMethod("DelayedArray", "repoArraySeed",
    function(seed) DelayedArray:::new_DelayedArray(seed, Class="repoArray")
)

### Works directly on a HDF5ArraySeed object, in which case it must be called
### with a single argument.
repoArray <- function(repopath, name)
{
    seed <- repoArraySeed(repopath, name)
    DelayedArray(seed)
}
