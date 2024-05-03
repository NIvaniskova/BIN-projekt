int parsefile(const char* fname, double* data, int* par_inputs, int* par_outputs)
{
    unsigned int inputs = 0;
    unsigned int outputs = 0;
    FILE    *fh = NULL;
    char    ch;   
    int res, lineno;
    double d;
    int lastseppos = 0;
    int values = 0;

    #define flushbuf() for (unsigned int i=0; i < inputs+outputs; i++) { if (data) { data[datarows] = temp[i];} /*printf("%d:%08X\n", datarows, temp[i]);*/ datarows++; }
    if (!(fh = fopen(fname,"r"))) {fprintf(stderr, "File '%s' not found\n", fname); return -1; }

    do {
       res = fscanf(fh, "%lf", &d);
       if (res == 0) 
       {
          ch = getc(fh);
          if (ch == '#') {
             while ((ch != '\n') && ((ch = getc(fh)) != EOF));
          } else if (ch == ':') { 
             if (inputs == 0) {
                inputs = values;
//                printf("inputs=%d ", inputs);                     
             } else if (outputs == 0) {
                outputs = values - 2*inputs;
//                printf("outputs=%d ", outputs);                     
             } else if (values - lastseppos != (inputs + outputs)) {
                 fprintf(stderr,"line %d: Invalid number of items",lineno);
                 return -1;
             }
             lastseppos = values;
          }
          
          res = 1;

       } else if (res > 0) {
          if (data != 0)
              data[values] = d;
//          printf("%lf\n", d);
          values++;
       }
    } while (res > 0);
    fclose(fh);

    if (par_inputs != 0) *par_inputs = inputs;
    if (par_outputs != 0) *par_outputs = outputs;

//    printf("\n");
    return values;
}
