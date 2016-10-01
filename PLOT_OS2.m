function PLOT_OS2(plot_type,ZAXIS,RAXIS,data,cmap_mat,fig_num,varargin)

if length(varargin)==2
    if strcmp(varargin{1},'cax')
        cax_val = varargin{2};
    end
else
    cax_val = 1;
end


figure(fig_num);

if strcmp(plot_type,'d1')
    plot(RAXIS,data,'k','linewidth',3); axis tight;
    xlabel('R [\mum]','fontsize',16);
    ylabel('n [1E15]','fontsize',16);
    title('Density Lineout','fontsize',16);
    set(gca,'fontsize',14);
end

if strcmp(plot_type,'density')
    imagesc(ZAXIS,RAXIS,data);
    axis xy;
    axis image;
    colormap(cmap_mat);
    caxis([-cax_val cax_val]);
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
    caxis([-emax*cax_val emax*cax_val]);
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
    caxis([-bmax*cax_val bmax*cax_val]);
    xlabel('Z [\mum]','fontsize',16);
    ylabel('R [\mum]','fontsize',16);
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String', 'E_r - B_{\theta} (MT/m)','fontsize',16);
    title('Focusing Field','fontsize',16);
    set(gca,'fontsize',14);
end