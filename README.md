# Transmission ratio distortion in species-wide crosses of ye

In this project, I am taking information that was available to me about *S. cerevisiae* yeast, based on genotypes and other data, to select which yeasts I want to cross <00_SelectCrosses>. 

Then, once the yeasts are crossed in the lab (not shown here), I get pooled sequencing reads of the segregants of the two crossed yeasts (per cross) and map them for further analysis in <01_Mapping>.

In <02_TRD> I take the mapped reads and extract where in the genome we see the allele frequencies deviate from the 0.5 (50%) expectation of Mendelian segregation, i.e. transmission ratio distortion.

Finally, for now, in <03_GenomicSignals>, I am currently investigating all sorts of different traits of the regions that show TRD signals. For this, I am leveraging the rich data sets of the Joseph Schacherer lab.

## Other folders

As I often do it, <Archive> contains notebooks etc that are not used anymore, either because they were done better or led nowhere etc. The files were `git mv`ed, so history can be kept.
    
In <RStudioScripts> I keep scripts that I used locally, rather than on a cluster. In this case, mostly wet lab related.
    
<scripts> contains R scripts that are used elsewhere, but is not added to anymore as I now keep scripts with their accompanying notebooks.
    
In <Shiny>, I am working on a shiny app to browse the crosses. Some data is implemented and more will come™️.