graph
[
node
[
	id 0
	v 1
	u1 1
	u2 1
]
node
[
	id 1
	v 2
	u1 2
	u2 2
]
node
[
	id 2
	v {2,3} 
	u1 3
	u2 3
]
node
[
	id 3
	v {1,2,3} 
	u1 1
	u2 {1,3} 
]
node
[
	id 4
	v {1,3} 
	u1 1
	u2 1
]
edge
[
	source 0
	target 1
]
edge
[
	source 0
	target 4
]
edge
[
	source 1
	target 2
]
edge
[
	source 2
	target 3
]
edge
[
	source 3
	target 4
]
]
