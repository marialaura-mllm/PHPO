functions{

  vector loglik_weibull(vector H0t, vector ht, vector status){

    vector[num_elements(status)] loglik;

    for(i in 1:num_elements(status)){loglik[i]=(status[i]*log(ht[i]))-H0t[i];}

    return loglik;
  }

  vector loglik_weibull_PO(vector H0t, vector h0t, vector status, vector lambda_PO){

    vector[num_elements(status)] loglik;

    for(i in 1:num_elements(status)){loglik[i]= status[i]*log((lambda_PO[i]*h0t[i])/((exp(-H0t[i]))+lambda_PO[i]*(1-exp(-H0t[i]))))-log(1+lambda_PO[i]*((1-exp(-H0t[i]))/(exp(-H0t[i]))));}

    return loglik;
  }

}
