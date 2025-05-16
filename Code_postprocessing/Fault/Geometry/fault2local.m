function vectors_local = fault2local(vectors_fault,strike,dip)

% fault2local.m takes as input n_vectors vectors in the Okada reference
% frame, where x is along strike, y is along dip (updip), and z is
% orthogonal to the fault plane, and transforms them in the local reference
% frame where x is along East, y along North, and z is orthogonal to the
% ground surface pointing outward the center of the Earth.
% -------------------------------------------------------------------------
% INPUT
% vectors_fault: strike, updip, and tensile components of the vectors to be
%                transformed; size: 3 x n_vectors.
% strike: strike angle of the fault, defined starting from the North
%         direction, increasing clockwise, until the strike direction is
%         reached; size: scalar.
% dip: dip angle of the fault, defined as the angle between the plane
%      parallel to the ground and the updip direction of the fault; size:
%      scalar.
% -------------------------------------------------------------------------
% OUTPUT
% vectors_local: East, North, and Vertical components of the transformed
%                vectors; size: 3 x n_vectors.
% -------------------------------------------------------------------------
% 
% Adriano Gualandi - 25 Aug 2016
% California Institute of Technology
% Geological and Planetary Science Division

% Let us call [x;y;z] the system from which we start, i.e. the fault
% reference system where x is along the strike direction, y is along the
% updip direction, and z is along the tensile direction.
% Let us call [x'';y'';z''] the system to which we end, i.e. the local
% geographic reference system where x'' is along the East direction, y'' is
% along the North direction, and z'' is along the Vertical outward
% direction.
% The transformation consists of two rotations. After the first rotation we
% will be in the reference system defined by the intermediate coordinates
% [x';y';z'].

% The x and x'' directions, i.e. the strike and East directions, lie on a
% plane parallel to the ground surface, i.e. parallel to the plane z''=0.
% For this reason as first rotation we align the z axis with the z'' axis,
% so that then we can rotate around z'(=z'') to align x' with x''. To do
% that we have to perform the rotation around the x axis, so that x'=x. The
% dip angle is defined as that angle that is formed between the plane
% parallel to the ground and the fault plane. It follows that an analogue
% definition can be the angle that goes from the vertical direction to the
% tensile direction. Thus, we have to perform a rotation of -dip to align
% the tensile axis z with the vertical axis z''. This rotation is performed
% around the strike direction x.
angle1 = -dip;

% After the first rotation we now have to align the x'=x strike direction
% with the East direction x''. We already know that now z has been aligned
% with the Vertical direction z'', so we do not have to modify z'. Thus, we
% can rotate around the axis z'=z''. The rotation must be of an angle equal
% to the angle that separates the strike and East directions. The angle
% between the strike direction and the East is given by (strike-90).
% Indeed, the strike angle is increasing clockwise from the North. If we
% subtract 90 degrees we get the angle between the East direction and the
% strike direction clockwise increasing, i.e. moving the East axis towards
% the strike direction. The rotation must be from the strike direction and
% the East direction of the same angle, that will now increase
% counterclockwise.
angle2 = (strike-90);

% The transformation matrices are calculated using the
% rotation_reference_system function.
% Rotation around the strike direction (x axis) of an angle equal to angle1
R1 = rotation_reference_system('x',angle1);
% Rotation around the Vertical direction (z'' axis) of an angle equal to
% angle2
R2 = rotation_reference_system('z',angle2);

% The overall transformation is:
% vectors_local = R2*R1*vectors_fault
% where vectors_local and vectors_fault are the vectors in the local and
% fault reference systems.
vectors_local = R2*R1*vectors_fault;