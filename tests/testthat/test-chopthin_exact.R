library("chopthin")
set.seed(1237123)


checkbasics <- function(res,N){
    expect_true(length(res$indices)==N);
    expect_true(length(res$weights)==N);
    expect_true(all(res$indices>0));
    expect_true(all(res$weights>0));
    expect_true(abs(sum(res$weights)-N)<1e-10); #weights normalised to N
}

test_that("Individual small scale tests",{
    expect_error(chopthin(rep(0,5),5));
    checkbasics(chopthin(1,N=1),1)
    
    for (i in c(100,101,105,150,197,198,199,200)){
        res <- chopthin(rep(1,100),N=i);
        checkbasics(res,i)
     }
})


test_that("samples",{
    set.seed(1230);
    w <- rexp(10);
    res <- replicate(1e4,chopthin(w,N=1)$indices)
    expect_true(chisq.test(table(res),p=w/sum(w))$p.value>1e-1)
    
    skip_on_cran()
    for (i in 1:10){
        w <- rexp(20);
        res <- replicate(1e5,chopthin(w,N=1)$indices)
        expect_true(chisq.test(table(res),p=w/sum(w))$p.value>1e-2)
    }
})
