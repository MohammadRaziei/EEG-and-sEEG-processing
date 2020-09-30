function plot_haedmodel_elec_topo(headmodel_brain, elec, potentials)
hold on; ft_plot_topo3d(elec.elecpos, potentials,'facealpha',0.5)
hold on; ft_plot_sens(elec,'style','blue');
hold on; ft_plot_mesh(headmodel_brain,'edgecolor','none','facecolor',[0.4 0.6 0.4],'facealpha',0.1);
