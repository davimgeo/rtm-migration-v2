
/*
void rtm()
{
  for (int isrc = 0; isrc < nsrc; isrc++) 
  {
    reset_fields();
    
    define_source_coordinates(isrc);

    for (int t = 1; t < nt-1; t++) 
    {
      forward_propagation(t);

      get_src_snaps(t);
    }

    for (int t = nt-1; t < tstop; t++) 
    {
      backward_propagation(t);

      accumulate_cross_correlation(t);
    }

   image_condition();

   // show_modeling_status
  } 
}
*/
