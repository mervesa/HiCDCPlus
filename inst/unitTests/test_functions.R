test_functions <- function() {
    checkEquals(get.chrs("Hsapiens","hg19"), paste0("chr",c(1:22,"X")))
}
