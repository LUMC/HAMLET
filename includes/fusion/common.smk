def get_fusioncatcher_outputs():
    # If no fusioncatcher data is specified, we do not run it and there are no
    # outputs
    if not settings["fusioncatcher_data"]:
        return list()

    fs = expand("{sample}/fusion/{sample}.fusioncatcher", sample=config["samples"])
    fuma = expand("{sample}/fusion/{sample}.fuma", sample=config["samples"])
    is_svg = expand("{sample}/fusion/{sample}.sf-isect-circos/fsnviz.svg", sample=config["samples"])
    is_png = expand("{sample}/fusion/{sample}.sf-isect-circos/fsnviz.png", sample=config["samples"])
    fc_svg = expand("{sample}/fusion/{sample}.fusioncatcher-circos/fsnviz.svg", sample=config["samples"])
    fc_png = expand("{sample}/fusion/{sample}.fusioncatcher-circos/fsnviz.png", sample=config["samples"])
    combined_svg = expand("{sample}/fusion/{sample}.fusions-combined.svg", sample=config["samples"])
    subset = expand("{sample}/fusion/{sample}.sf-isect", sample=config["samples"])
    return fs + fuma + is_svg + is_png + fc_svg + fc_png + combined_svg + subset
