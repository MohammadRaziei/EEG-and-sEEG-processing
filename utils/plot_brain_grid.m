function plot_brain_grid(leadfield)
resolution = max(abs(diff(leadfield.pos(1:2,:))));

plot3(leadfield.pos(:, 1), ...
  leadfield.pos(:, 2), ...
  leadfield.pos(:, 3), 'o', ...
  'markersize', 2*resolution, 'Color', '#555',...
  'markeredgecolor','#999')
hold on
plot3(leadfield.pos(leadfield.inside, 1), ...
  leadfield.pos(leadfield.inside, 2), ...
  leadfield.pos(leadfield.inside, 3), 'o', ...
  'markersize', 2*resolution, 'markerfacecolor', '#eee',...
  'markeredgecolor','#999')


end


