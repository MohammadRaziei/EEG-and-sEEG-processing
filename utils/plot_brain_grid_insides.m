function plot_brain_grid_insides(leadfield, resolution)
if nargin < 2 , resolution = max(abs(diff(leadfield.pos(1:2,:)))); end

plot3(leadfield.pos(leadfield.inside, 1), ...
    leadfield.pos(leadfield.inside, 2), ...
    leadfield.pos(leadfield.inside, 3), 'o', ...
    'markersize', 2*resolution,...
    'markeredgecolor','#999')
end

