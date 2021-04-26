import torch
import torchvision

print(torch.__version__)
a = [1.0, 2.0, 3.0]
print(a[2])
a = torch.ones(3)
print(a)
float(a[1])
a[1]
points = torch.zeros(6)
print(points)
points = torch.tensor([1.0, 4.0, 2.0, 1.0, 3.0, 5.0])
print(points)
points = torch.tensor([[1.0, 4.0], [2.0, 1.0], [3.0, 5.0]])
print(points)
points.shape
torch.__version__
points = torch.zeros(3,2)
print(points)
points = torch.FloatTensor([[1.0, 4.0], [2.0, 1.0], [3.0, 5.0]])
print(points)
print(points[2])
points.storage()[5]
point2 = points[1]
point2
point2.storage_offset()
point2.size()
point1 = points[0]
point1.storage_offset()
point1.size()
point3 = points[2]
point3
point3.storage_offset()
point3.size()
point3.stride()
point2.stride()
points.stride()
points.storage_offset()
point2 = points[1].clone()
point2[0] = 10.0
points
point2 = points[1]
point2[0] = 10.0
points
points_t = points.t()
points_t
points_t[1,1] = 2
points.size()
points.stride()
points_t.size()
points_t.stride()

some = torch.ones(3, 4, 5)
some[1, 2, 4]
some[0, 2]
some.stride()
some_t = some.transpose(2,1)
some_t
some
some_t.storage_offset()
point.storage()
points = torch.tensor([[1.0, 4.0], [2.0, 1.0], [3.0, 5.0]])
points_t = points.t()
points_t
points_t.storage()
points_t.stride()
points.stride()
points_t_cont = points_t.contiguous()
points_t_cont.stride()
points_t_cont
points.stride()
points
points_t_cont[1, 1]=10.0
points = torch.randn(10, 2)
points.type()
short_points = points.type(torch.short)
short_points
short_points[1, 1] = 1
short_points
points
some_list = list(range(6))
some_list
points
points[1:, 0]
