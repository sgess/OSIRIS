function PLOT_OS2(plot_type,ZAXIS,RAXIS,data,cmap_mat,fig_num)

figure(fig_num);

if strcmp(plot_type,'density')
    imagesc(ZAXIS,RAXIS,data);
    axis xy;
    axis image;
    colormap(cmap_mat);
    caxis([-1 1]);
    xlabel('Z [\mum]','fontsize',16);
    ylabel('R [\mum]','fontsize',16);
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String', 'n_0','fontsize',16);
    title('Density','fontsize',16);
    set(gca,'fontsize',14);
end

if strcmp(plot_type,'ez2')
    emax = max(abs(data(:)));
    imagesc(ZAXIS,RAXIS,data);
    axis xy;
    axis image;
    colormap(cmap_mat);
    caxis([-emax emax]);
    xlabel('Z [\mum]','fontsize',16);
    ylabel('R [\mum]','fontsize',16);
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String', 'E_z (GV/m)','fontsize',16);
    title('Longitudinal Field','fontsize',16);
    set(gca,'fontsize',14);
end

if strcmp(plot_type,'ez1')
    plot(ZAXIS,data,'k','linewidth',3); axis tight;
    xlabel('Z [\mum]','fontsize',16);
    ylabel('E_z (GV/m)','fontsize',16);
    title('Longitudinal Field','fontsize',16);
    set(gca,'fontsize',14);
end

if strcmp(plot_type,'fr2')
    bmax = max(abs(data(:)));
    imagesc(ZAXIS,RAXIS,data);
    axis xy;
    axis image;
    colormap(cmap_mat);
    caxis([-bmax bmax]);
    xlabel('Z [\mum]','fontsize',16);
    ylabel('R [\mum]','fontsize',16);
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String', 'E_r - B_{\theta} (MT/m)','fontsize',16);
    title('Focusing Field','fontsize',16);
    set(gca,'fontsize',14);
end