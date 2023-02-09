stringr_count_bases = function(string_column,ref_column)
    {
    #if(nchar(string) > 0 & nchar(ref) > 0)
     #   {
        

        a=str_count(string_column, fixed("a"))
        a=a+str_count(string_column, fixed("A"))

        t=str_count(string_column, fixed("t"))
        t=t+str_count(string_column, fixed("T"))

        g=str_count(string_column, fixed("g"))
        g=g+str_count(string_column, fixed("G"))

        c=str_count(string_column, fixed("c"))
        c=c+str_count(string_column, fixed("C"))
    
        ref=str_count(string_column, fixed(","))
        ref=ref+str_count(string_column, fixed("."))

        a = a + (ref * (ref_column == 'a'))
        a = a + (ref * (ref_column == 'A'))
    
        t = t + (ref * (ref_column == 't'))
        t = t + (ref * (ref_column == 'T'))
    
        g = g + (ref * (ref_column == 'g'))
        g = g + (ref * (ref_column == 'G'))
    
        c = c + (ref * (ref_column == 'c'))
        c = c + (ref * (ref_column == 'C'))


        return(cbind(t,a,g,c))
      #  }
   # else
     #   {
       # return(cbind(NA,NA,NA,NA))
  #  }
    
}

ascii_score_to_numeric = function(ascii)
    {
    return(mean(utf8ToInt(ascii)-33))
}