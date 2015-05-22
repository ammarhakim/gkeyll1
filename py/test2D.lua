-- Input file for Poisson bracket operator

-- grid on which equations are to be solved
grid = Grid.RectCart2D {
   lower = {0, 0},
   upper = {1, 2},
   cells = {8, 16},
}

function pulse(x,y,z)
   local xc, yc = 0.5, 0.5
   local r2 = (x-xc)^2 + (y-yc)^2
   return math.exp(-75*r2)
end

basis1 = NodalFiniteElement2D.LagrangeTensor {
   onGrid = grid,
   polyOrder = 1,
   nodeLocation = "lobatto",
}
q1 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis1:numNodes(),
   ghost = {2, 2},
}
initField1 = Updater.ProjectOnNodalBasis2D {
   onGrid = grid,
   basis = basis1,
   shareCommonNodes = false,
   evaluate = function (x,y,z,t)
		 return pulse(x,y,z)
	      end
}
initField1:setOut( {q1} )
initField1:advance(0.0)
q1:write("lag-2D-p1.h5")

basis2 = NodalFiniteElement2D.LagrangeTensor {
   onGrid = grid,
   polyOrder = 2,
   nodeLocation = "lobatto",
}
q2 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis2:numNodes(),
   ghost = {2, 2},
}
initField2 = Updater.ProjectOnNodalBasis2D {
   onGrid = grid,
   basis = basis2,
   shareCommonNodes = false,
   evaluate = function (x,y,z,t)
		 return pulse(x,y,z)
	      end
}
initField2:setOut( {q2} )
initField2:advance(0.0)
q2:write("lag-2D-p2.h5")

basis3 = NodalFiniteElement2D.LagrangeTensor {
   onGrid = grid,
   polyOrder = 3,
   nodeLocation = "lobatto",
}
q3 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis3:numNodes(),
   ghost = {2, 2},
}
initField3 = Updater.ProjectOnNodalBasis2D {
   onGrid = grid,
   basis = basis3,
   shareCommonNodes = false,
   evaluate = function (x,y,z,t)
		 return pulse(x,y,z)
	      end
}
initField3:setOut( {q3} )
initField3:advance(0.0)
q3:write("lag-2D-p3.h5")

basis4 = NodalFiniteElement2D.LagrangeTensor {
   onGrid = grid,
   polyOrder = 4,
   nodeLocation = "lobatto",
}
q4 = DataStruct.Field2D {
   onGrid = grid,
   numComponents = basis4:numNodes(),
   ghost = {2, 2},
}
initField4 = Updater.ProjectOnNodalBasis2D {
   onGrid = grid,
   basis = basis4,
   shareCommonNodes = false,
   evaluate = function (x,y,z,t)
		 return pulse(x,y,z)
	      end
}
initField4:setOut( {q4} )
initField4:advance(0.0)
q4:write("lag-2D-p4.h5")