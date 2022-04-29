functions {
  vector project_smplx(matrix x, int k, int t) {
    row_vector[rows(x)] y = rep_row_vector(0, rows(x));
    y[k] = 1;
    for(j in 1:t) {
      y = y * x;
    }
    return to_vector(y);
  }

}


model {
  
  
}
