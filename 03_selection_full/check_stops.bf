RequireVersion ("2.5.50");
fprintf (stdout, "Script starting...\n");

base_path = "/Users/sergei/Dropbox/Work/Collaborations/Y-chromosome/2026/round6_export_20260126/";

file_list = {};
file_list [0] = "BPY2/BPY2.unique.nexus";
file_list [1] = "CDY/CDY.unique.nexus";
file_list [2] = "CDY_CDYL/CDY_CDYL.unique.nexus";
file_list [3] = "DAZ_DAZL_repeats/DAZ_DAZL_repeats.unique.nexus";
file_list [4] = "DAZ_DAZL_RRM/DAZ_DAZL_RRM.unique.nexus";
file_list [5] = "DAZ_repeats/DAZ_repeats.unique.nexus";
file_list [6] = "DAZ_RRM/DAZ_RRM.unique.nexus";
file_list [7] = "HSFY/HSFY.unique.nexus";
file_list [8] = "RBMY/RBMY.unique.nexus";
file_list [9] = "RBMY_RBMX/RBMY_RBMX.unique.nexus";
file_list [10] = "TSPY/TSPY.unique.nexus";
file_list [11] = "VCY/VCY.unique.nexus";

for (file_id = 0; file_id < Abs (file_list); file_id += 1) {
    file_path = base_path + file_list[file_id];
    fprintf (stdout, "Checking ", file_path, "...\n");
    
    DataSet ds = ReadDataFile (file_path);
    
    if (ds.sites % 3 != 0) {
        fprintf (stdout, "  [WARNING] Alignment length (", ds.sites, ") is not a multiple of 3.\n");
    }
    
    DataSetFilter codon_data = CreateFilter (ds, 3, "", "", "TAA,TAG,TGA");
    GetInformation (seq_names, codon_data);
    
    for (i = 0; i < codon_data.species; i += 1) {
        GetString (raw_seq, ds, i);
        stop_count = 0;
        len = Abs(raw_seq);
        
        for (k = 0; k < len - 3; k += 3) { 
            codon = raw_seq[k] + raw_seq[k+1] + raw_seq[k+2];
            // Simple case-insensitive check
            if (codon == "TAA" || codon == "TAG" || codon == "TGA" || 
                codon == "taa" || codon == "tag" || codon == "tga") {
                stop_count += 1;
            }
        }
        
        if (stop_count > 0) {
             fprintf (stdout, "  [WARNING] Sequence '", seq_names[i], "' has ", stop_count, " internal premature stop codons.\n");
        }
    }
}