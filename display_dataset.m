function r = display_dataset(dataset)
data = sum(permute(dataset,[3,1,2]))
image(data);