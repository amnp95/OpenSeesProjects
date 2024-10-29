

def infowriter(info,meshdir):
    """
    Write the information dictionary to a tcl file.

    Parameters
    ----------
    info : dict
        Information dictionary.
    meshdir : str
        Directory where the tcl file will be written
            
    Returns
    -------
    None

    """
    # ============================================================================
    # Create a tcl file
    # ============================================================================
    f = open(f"{meshdir}/Modelinfo.tcl", "w")
    f.write(f"wipe\n")
    f.write(f"# ============================================================================\n")
    f.write(f"# Cores Information\n")
    f.write(f"# ============================================================================\n")
    f.write(f"set regcores       {info['regcores']}\n")
    f.write(f"set pmlcores       {info['pmlcores']}\n")
    f.write(f"set drmcores       {info['drmcores']}\n")
    f.write(f"set structurecores {info['structurecores']}\n")
    f.write(f"set AnalysisType   {info['AnalysisType']}\n")


    f.write(f"# ============================================================================\n")
    f.write(f"# Model Information\n")
    f.write(f"# ============================================================================\n")
    f.write(f"set StructureType     \"{info['StructureType']}\"\n")
    f.write(f"set NStory            {info['NStory']}\n")
    f.write(f"set NBay              {info['NBay']}\n")
    f.write(f"set NBayZ             {info['NBayZ']}\n")
    f.write(f"set StartNodeX       {info['StartNodeX']}\n")
    f.write(f"set StartNodeY       {info['StartNodeY']}\n")
    f.write(f"set StartNodeZ       {info['StartNodeZ']}\n")
    f.write(f"set LCol              {info['LCol']}\n")
    f.write(f"set LBeam             {info['LBeam']}\n")
    f.write(f"set LGird             {info['LGird']}\n")
    f.write(f"set SectionType       {info['SectionType']}\n")
    f.write(f"set HaveStructure     \"{info['HaveStructure']}\"\n")
    f.write(f"# ============================================================================\n")
    f.write(f"# Soil Information\n")
    f.write(f"# ============================================================================\n")
    f.write(f"set dx                {info['dx']}\n")
    f.write(f"set dy                {info['dy']}\n")
    f.write(f"set dz                {info['dz']}\n")
    f.write(f"set llx               {info['llx']}\n")
    f.write(f"set lly               {info['lly']}\n")
    f.write(f"set llz               {info['llz']}\n")
    f.write(f"set drmthicknessx     {info['drmthicknessx']}\n")
    f.write(f"set drmthicknessy     {info['drmthicknessy']}\n")
    f.write(f"set drmthicknessz     {info['drmthicknessz']}\n")
    f.write(f"set numdrmlayers      {info['numdrmlayers']}\n")
    f.write(f"set lx                {info['lx']}\n")
    f.write(f"set ly                {info['ly']}\n")
    f.write(f"set lz                {info['lz']}\n")
    f.write(f"set nx                {info['nx']}\n")
    f.write(f"set ny                {info['ny']}\n")
    f.write(f"set nz                {info['nz']}\n")
    f.write(f"# ============================================================================\n")
    f.write(f"# PML information\n")
    f.write(f"# ============================================================================\n")
    f.write(f"set AbsorbingElements \"{info['AbsorbingElements']}\"\n")
    f.write(f"set numpmllayers      {info['numpmllayers']}\n")
    f.write(f"set pmlthicknessx     {info['pmlthicknessx']}\n")
    f.write(f"set pmlthicknessy     {info['pmlthicknessy']}\n")
    f.write(f"set pmlthicknessz     {info['pmlthicknessz']}\n")
    f.write(f"set pmltotalthickness {info['pmltotalthickness']}\n")
    f.write(f"set HaveAbsorbingElements \"{info['HaveAbsorbingElements']}\"\n")
    f.write(f"# ============================================================================\n")
    f.write(f"# General information\n")
    f.write(f"# ============================================================================\n")
    f.write(f"set meshdir           \"{info['meshdir']}\"\n")
    f.write(f"set outputdir         \"{info['outputdir']}\"\n")
    f.write(f"set DRMFile           \"{info['DRMFile']}\"\n")
    f.write(f"# ============================================================================\n")
    f.write(f"# Embedding foundation\n")
    f.write(f"# ============================================================================\n")
    f.write(f"set EmbeddingFoundation \"{info['EmbeddingFoundation']}\"\n")
    info2 = info['EmbeddedFoundation']
    f.write(f"set EmbeddedFoundation [dict create xmax {info2['xmax']} xmin {info2['xmin']} ymax {info2['ymax']} ymin {info2['ymin']} zmax {info2['zmax']} zmin {info2['zmin']}]\n")
    f.write(f"# ============================================================================\n")
    f.write(f"# Fondation information\n")
    f.write(f"# ============================================================================\n")
    f.write(f"set HaveFoundation \"{info['HaveFoundation']}\"\n")
    f.write("set foundationBlocks {}\n")
    for i,block in enumerate(info['foundationBlocks']):
        f.write(f"lappend foundationBlocks [dict create matTag {block['matTag']} xmax {block['xmax']} xmin {block['xmin']} ymax {block['ymax']} ymin {block['ymin']} zmax {block['zmax']} zmin {block['zmin']} Xmeshsize {block['Xmeshsize']} Ymeshsize {block['Ymeshsize']} Zmeshsize {block['Zmeshsize']}]\n")
    f.write(f"# ============================================================================\n")
    f.write(f"# Piles information\n")
    f.write(f"# ============================================================================\n")
    f.write(f"set HavePiles \"{info['HavePiles']}\"\n")
    f.write("set pilelist {}\n")
    for i,pile in enumerate(info['pilelist']):
        f.write(f"lappend pilelist [dict create xtop {pile['xtop']} ytop {pile['ytop']} ztop {pile['ztop']} xbottom {pile['xbottom']} ybottom {pile['ybottom']} zbottom {pile['zbottom']} numberofElements {pile['numberofElements']}]\n")
    f.write(f"# ============================================================================\n")
    f.write(f"# cells and nodes information\n")
    f.write(f"# ============================================================================\n")
    f.write(f"set soilfoundation_num_cells {info['soilfoundation_num_cells']}\n")
    f.write(f"set soilfoundation_num_nodes {info['soilfoundation_num_points']}\n")
    f.close()

    return None






