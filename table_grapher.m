function table_grapher(table_name) 
for i=1:size(table_name.runtime)
    loglog(table_name.error(i),table_name.runtime(i),table_name.markers{i})
    hold on;
end
hold off;
xlabel 'reconstruction error'
ylabel 'runtime'
%lh = legend(table_name.Properties.RowNames);
%set(lh,'location','northeastoutside');
axis([0.001,0.1,0.01,10000])

end