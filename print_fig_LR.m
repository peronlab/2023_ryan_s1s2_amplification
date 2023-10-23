function print_fig_LR(fh, output_filename)
    print_on = 1;
 
    if (~print_on) ; return ; end
     
    print_dir = '~/Desktop/Paper_figs/fig5';
    if (~exist(print_dir, 'dir'))
        mkdir(print_dir);
    end
 
    set(fh,'PaperPositionMode','auto');
 
    ppos = get(fh,'PaperPosition');
    su = get(fh,'Units');
    pu = get(fh,'PaperUnits'); 
    set(fh,'Units',pu);
    spos = get(fh,'Position');
    set(fh,'Position',[spos(1) spos(2) ppos(3) ppos(4)])
    set(fh,'Units',su)
    set(fh,'renderer','Painters');
 
    % for Illustrator
    output_full_path = [print_dir filesep output_filename '.eps'];
    print (fh, output_full_path ,'-dpsc2', '-painters', '-noui');
 
    % for Inkscape
%     output_full_path = [print_dir filesep output_filename '.svg'];
%     print (fh, output_full_path ,'-dsvg', '-painters', '-noui');