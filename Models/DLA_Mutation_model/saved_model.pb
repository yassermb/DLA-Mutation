??6
??
8
Const
output"dtype"
valuetensor"
dtypetype

NoOp
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
@
ReadVariableOp
resource
value"dtype"
dtypetype?
?
StatefulPartitionedCall
args2Tin
output2Tout"
Tin
list(type)("
Tout
list(type)("	
ffunc"
configstring "
config_protostring "
executor_typestring ?
?
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape"#
allowed_deviceslist(string)
 ?"serve*2.3.02v2.3.0-2-gee598066c48??.
x
dense_2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*
shared_namedense_2/kernel
q
"dense_2/kernel/Read/ReadVariableOpReadVariableOpdense_2/kernel*
_output_shapes

:*
dtype0
p
dense_2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_namedense_2/bias
i
 dense_2/bias/Read/ReadVariableOpReadVariableOpdense_2/bias*
_output_shapes
:*
dtype0
f
	Adam/iterVarHandleOp*
_output_shapes
: *
dtype0	*
shape: *
shared_name	Adam/iter
_
Adam/iter/Read/ReadVariableOpReadVariableOp	Adam/iter*
_output_shapes
: *
dtype0	
j
Adam/beta_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/beta_1
c
Adam/beta_1/Read/ReadVariableOpReadVariableOpAdam/beta_1*
_output_shapes
: *
dtype0
j
Adam/beta_2VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/beta_2
c
Adam/beta_2/Read/ReadVariableOpReadVariableOpAdam/beta_2*
_output_shapes
: *
dtype0
h

Adam/decayVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name
Adam/decay
a
Adam/decay/Read/ReadVariableOpReadVariableOp
Adam/decay*
_output_shapes
: *
dtype0
x
Adam/learning_rateVarHandleOp*
_output_shapes
: *
dtype0*
shape: *#
shared_nameAdam/learning_rate
q
&Adam/learning_rate/Read/ReadVariableOpReadVariableOpAdam/learning_rate*
_output_shapes
: *
dtype0
?
conv3d/kernelVarHandleOp*
_output_shapes
: *
dtype0* 
shape:?*
shared_nameconv3d/kernel
|
!conv3d/kernel/Read/ReadVariableOpReadVariableOpconv3d/kernel*+
_output_shapes
:?*
dtype0
n
conv3d/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameconv3d/bias
g
conv3d/bias/Read/ReadVariableOpReadVariableOpconv3d/bias*
_output_shapes
:*
dtype0
?
batch_normalization/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:**
shared_namebatch_normalization/gamma
?
-batch_normalization/gamma/Read/ReadVariableOpReadVariableOpbatch_normalization/gamma*
_output_shapes
:*
dtype0
?
batch_normalization/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*)
shared_namebatch_normalization/beta
?
,batch_normalization/beta/Read/ReadVariableOpReadVariableOpbatch_normalization/beta*
_output_shapes
:*
dtype0
?
conv3d_1/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:* 
shared_nameconv3d_1/kernel

#conv3d_1/kernel/Read/ReadVariableOpReadVariableOpconv3d_1/kernel**
_output_shapes
:*
dtype0
r
conv3d_1/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameconv3d_1/bias
k
!conv3d_1/bias/Read/ReadVariableOpReadVariableOpconv3d_1/bias*
_output_shapes
:*
dtype0
?
batch_normalization_1/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*,
shared_namebatch_normalization_1/gamma
?
/batch_normalization_1/gamma/Read/ReadVariableOpReadVariableOpbatch_normalization_1/gamma*
_output_shapes
:*
dtype0
?
batch_normalization_1/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*+
shared_namebatch_normalization_1/beta
?
.batch_normalization_1/beta/Read/ReadVariableOpReadVariableOpbatch_normalization_1/beta*
_output_shapes
:*
dtype0
?
conv3d_2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:* 
shared_nameconv3d_2/kernel

#conv3d_2/kernel/Read/ReadVariableOpReadVariableOpconv3d_2/kernel**
_output_shapes
:*
dtype0
r
conv3d_2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameconv3d_2/bias
k
!conv3d_2/bias/Read/ReadVariableOpReadVariableOpconv3d_2/bias*
_output_shapes
:*
dtype0
?
batch_normalization_2/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*,
shared_namebatch_normalization_2/gamma
?
/batch_normalization_2/gamma/Read/ReadVariableOpReadVariableOpbatch_normalization_2/gamma*
_output_shapes
:*
dtype0
?
batch_normalization_2/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*+
shared_namebatch_normalization_2/beta
?
.batch_normalization_2/beta/Read/ReadVariableOpReadVariableOpbatch_normalization_2/beta*
_output_shapes
:*
dtype0
?
conv3d_3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:* 
shared_nameconv3d_3/kernel

#conv3d_3/kernel/Read/ReadVariableOpReadVariableOpconv3d_3/kernel**
_output_shapes
:*
dtype0
r
conv3d_3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameconv3d_3/bias
k
!conv3d_3/bias/Read/ReadVariableOpReadVariableOpconv3d_3/bias*
_output_shapes
:*
dtype0
?
batch_normalization_3/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*,
shared_namebatch_normalization_3/gamma
?
/batch_normalization_3/gamma/Read/ReadVariableOpReadVariableOpbatch_normalization_3/gamma*
_output_shapes
:*
dtype0
?
batch_normalization_3/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:*+
shared_namebatch_normalization_3/beta
?
.batch_normalization_3/beta/Read/ReadVariableOpReadVariableOpbatch_normalization_3/beta*
_output_shapes
:*
dtype0
x
layer1/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:
?
?*
shared_namelayer1/kernel
q
!layer1/kernel/Read/ReadVariableOpReadVariableOplayer1/kernel* 
_output_shapes
:
?
?*
dtype0
o
layer1/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:?*
shared_namelayer1/bias
h
layer1/bias/Read/ReadVariableOpReadVariableOplayer1/bias*
_output_shapes	
:?*
dtype0
u
dense/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	?d*
shared_namedense/kernel
n
 dense/kernel/Read/ReadVariableOpReadVariableOpdense/kernel*
_output_shapes
:	?d*
dtype0
l

dense/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:d*
shared_name
dense/bias
e
dense/bias/Read/ReadVariableOpReadVariableOp
dense/bias*
_output_shapes
:d*
dtype0
x
dense_1/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:d
*
shared_namedense_1/kernel
q
"dense_1/kernel/Read/ReadVariableOpReadVariableOpdense_1/kernel*
_output_shapes

:d
*
dtype0
p
dense_1/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*
shared_namedense_1/bias
i
 dense_1/bias/Read/ReadVariableOpReadVariableOpdense_1/bias*
_output_shapes
:
*
dtype0
?
batch_normalization/moving_meanVarHandleOp*
_output_shapes
: *
dtype0*
shape:*0
shared_name!batch_normalization/moving_mean
?
3batch_normalization/moving_mean/Read/ReadVariableOpReadVariableOpbatch_normalization/moving_mean*
_output_shapes
:*
dtype0
?
#batch_normalization/moving_varianceVarHandleOp*
_output_shapes
: *
dtype0*
shape:*4
shared_name%#batch_normalization/moving_variance
?
7batch_normalization/moving_variance/Read/ReadVariableOpReadVariableOp#batch_normalization/moving_variance*
_output_shapes
:*
dtype0
?
!batch_normalization_1/moving_meanVarHandleOp*
_output_shapes
: *
dtype0*
shape:*2
shared_name#!batch_normalization_1/moving_mean
?
5batch_normalization_1/moving_mean/Read/ReadVariableOpReadVariableOp!batch_normalization_1/moving_mean*
_output_shapes
:*
dtype0
?
%batch_normalization_1/moving_varianceVarHandleOp*
_output_shapes
: *
dtype0*
shape:*6
shared_name'%batch_normalization_1/moving_variance
?
9batch_normalization_1/moving_variance/Read/ReadVariableOpReadVariableOp%batch_normalization_1/moving_variance*
_output_shapes
:*
dtype0
?
!batch_normalization_2/moving_meanVarHandleOp*
_output_shapes
: *
dtype0*
shape:*2
shared_name#!batch_normalization_2/moving_mean
?
5batch_normalization_2/moving_mean/Read/ReadVariableOpReadVariableOp!batch_normalization_2/moving_mean*
_output_shapes
:*
dtype0
?
%batch_normalization_2/moving_varianceVarHandleOp*
_output_shapes
: *
dtype0*
shape:*6
shared_name'%batch_normalization_2/moving_variance
?
9batch_normalization_2/moving_variance/Read/ReadVariableOpReadVariableOp%batch_normalization_2/moving_variance*
_output_shapes
:*
dtype0
?
!batch_normalization_3/moving_meanVarHandleOp*
_output_shapes
: *
dtype0*
shape:*2
shared_name#!batch_normalization_3/moving_mean
?
5batch_normalization_3/moving_mean/Read/ReadVariableOpReadVariableOp!batch_normalization_3/moving_mean*
_output_shapes
:*
dtype0
?
%batch_normalization_3/moving_varianceVarHandleOp*
_output_shapes
: *
dtype0*
shape:*6
shared_name'%batch_normalization_3/moving_variance
?
9batch_normalization_3/moving_variance/Read/ReadVariableOpReadVariableOp%batch_normalization_3/moving_variance*
_output_shapes
:*
dtype0
^
totalVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nametotal
W
total/Read/ReadVariableOpReadVariableOptotal*
_output_shapes
: *
dtype0
^
countVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_namecount
W
count/Read/ReadVariableOpReadVariableOpcount*
_output_shapes
: *
dtype0
?
Adam/dense_2/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*&
shared_nameAdam/dense_2/kernel/m

)Adam/dense_2/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense_2/kernel/m*
_output_shapes

:*
dtype0
~
Adam/dense_2/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*$
shared_nameAdam/dense_2/bias/m
w
'Adam/dense_2/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense_2/bias/m*
_output_shapes
:*
dtype0
?
Adam/conv3d/kernel/mVarHandleOp*
_output_shapes
: *
dtype0* 
shape:?*%
shared_nameAdam/conv3d/kernel/m
?
(Adam/conv3d/kernel/m/Read/ReadVariableOpReadVariableOpAdam/conv3d/kernel/m*+
_output_shapes
:?*
dtype0
|
Adam/conv3d/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*#
shared_nameAdam/conv3d/bias/m
u
&Adam/conv3d/bias/m/Read/ReadVariableOpReadVariableOpAdam/conv3d/bias/m*
_output_shapes
:*
dtype0
?
 Adam/batch_normalization/gamma/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*1
shared_name" Adam/batch_normalization/gamma/m
?
4Adam/batch_normalization/gamma/m/Read/ReadVariableOpReadVariableOp Adam/batch_normalization/gamma/m*
_output_shapes
:*
dtype0
?
Adam/batch_normalization/beta/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*0
shared_name!Adam/batch_normalization/beta/m
?
3Adam/batch_normalization/beta/m/Read/ReadVariableOpReadVariableOpAdam/batch_normalization/beta/m*
_output_shapes
:*
dtype0
?
Adam/conv3d_1/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*'
shared_nameAdam/conv3d_1/kernel/m
?
*Adam/conv3d_1/kernel/m/Read/ReadVariableOpReadVariableOpAdam/conv3d_1/kernel/m**
_output_shapes
:*
dtype0
?
Adam/conv3d_1/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*%
shared_nameAdam/conv3d_1/bias/m
y
(Adam/conv3d_1/bias/m/Read/ReadVariableOpReadVariableOpAdam/conv3d_1/bias/m*
_output_shapes
:*
dtype0
?
"Adam/batch_normalization_1/gamma/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*3
shared_name$"Adam/batch_normalization_1/gamma/m
?
6Adam/batch_normalization_1/gamma/m/Read/ReadVariableOpReadVariableOp"Adam/batch_normalization_1/gamma/m*
_output_shapes
:*
dtype0
?
!Adam/batch_normalization_1/beta/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*2
shared_name#!Adam/batch_normalization_1/beta/m
?
5Adam/batch_normalization_1/beta/m/Read/ReadVariableOpReadVariableOp!Adam/batch_normalization_1/beta/m*
_output_shapes
:*
dtype0
?
Adam/conv3d_2/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*'
shared_nameAdam/conv3d_2/kernel/m
?
*Adam/conv3d_2/kernel/m/Read/ReadVariableOpReadVariableOpAdam/conv3d_2/kernel/m**
_output_shapes
:*
dtype0
?
Adam/conv3d_2/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*%
shared_nameAdam/conv3d_2/bias/m
y
(Adam/conv3d_2/bias/m/Read/ReadVariableOpReadVariableOpAdam/conv3d_2/bias/m*
_output_shapes
:*
dtype0
?
"Adam/batch_normalization_2/gamma/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*3
shared_name$"Adam/batch_normalization_2/gamma/m
?
6Adam/batch_normalization_2/gamma/m/Read/ReadVariableOpReadVariableOp"Adam/batch_normalization_2/gamma/m*
_output_shapes
:*
dtype0
?
!Adam/batch_normalization_2/beta/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*2
shared_name#!Adam/batch_normalization_2/beta/m
?
5Adam/batch_normalization_2/beta/m/Read/ReadVariableOpReadVariableOp!Adam/batch_normalization_2/beta/m*
_output_shapes
:*
dtype0
?
Adam/conv3d_3/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*'
shared_nameAdam/conv3d_3/kernel/m
?
*Adam/conv3d_3/kernel/m/Read/ReadVariableOpReadVariableOpAdam/conv3d_3/kernel/m**
_output_shapes
:*
dtype0
?
Adam/conv3d_3/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*%
shared_nameAdam/conv3d_3/bias/m
y
(Adam/conv3d_3/bias/m/Read/ReadVariableOpReadVariableOpAdam/conv3d_3/bias/m*
_output_shapes
:*
dtype0
?
"Adam/batch_normalization_3/gamma/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*3
shared_name$"Adam/batch_normalization_3/gamma/m
?
6Adam/batch_normalization_3/gamma/m/Read/ReadVariableOpReadVariableOp"Adam/batch_normalization_3/gamma/m*
_output_shapes
:*
dtype0
?
!Adam/batch_normalization_3/beta/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*2
shared_name#!Adam/batch_normalization_3/beta/m
?
5Adam/batch_normalization_3/beta/m/Read/ReadVariableOpReadVariableOp!Adam/batch_normalization_3/beta/m*
_output_shapes
:*
dtype0
?
Adam/layer1/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
?
?*%
shared_nameAdam/layer1/kernel/m

(Adam/layer1/kernel/m/Read/ReadVariableOpReadVariableOpAdam/layer1/kernel/m* 
_output_shapes
:
?
?*
dtype0
}
Adam/layer1/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:?*#
shared_nameAdam/layer1/bias/m
v
&Adam/layer1/bias/m/Read/ReadVariableOpReadVariableOpAdam/layer1/bias/m*
_output_shapes	
:?*
dtype0
?
Adam/dense/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:	?d*$
shared_nameAdam/dense/kernel/m
|
'Adam/dense/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense/kernel/m*
_output_shapes
:	?d*
dtype0
z
Adam/dense/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:d*"
shared_nameAdam/dense/bias/m
s
%Adam/dense/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense/bias/m*
_output_shapes
:d*
dtype0
?
Adam/dense_1/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:d
*&
shared_nameAdam/dense_1/kernel/m

)Adam/dense_1/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense_1/kernel/m*
_output_shapes

:d
*
dtype0
~
Adam/dense_1/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*$
shared_nameAdam/dense_1/bias/m
w
'Adam/dense_1/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense_1/bias/m*
_output_shapes
:
*
dtype0
?
Adam/dense_2/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*&
shared_nameAdam/dense_2/kernel/v

)Adam/dense_2/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense_2/kernel/v*
_output_shapes

:*
dtype0
~
Adam/dense_2/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*$
shared_nameAdam/dense_2/bias/v
w
'Adam/dense_2/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense_2/bias/v*
_output_shapes
:*
dtype0
?
Adam/conv3d/kernel/vVarHandleOp*
_output_shapes
: *
dtype0* 
shape:?*%
shared_nameAdam/conv3d/kernel/v
?
(Adam/conv3d/kernel/v/Read/ReadVariableOpReadVariableOpAdam/conv3d/kernel/v*+
_output_shapes
:?*
dtype0
|
Adam/conv3d/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*#
shared_nameAdam/conv3d/bias/v
u
&Adam/conv3d/bias/v/Read/ReadVariableOpReadVariableOpAdam/conv3d/bias/v*
_output_shapes
:*
dtype0
?
 Adam/batch_normalization/gamma/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*1
shared_name" Adam/batch_normalization/gamma/v
?
4Adam/batch_normalization/gamma/v/Read/ReadVariableOpReadVariableOp Adam/batch_normalization/gamma/v*
_output_shapes
:*
dtype0
?
Adam/batch_normalization/beta/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*0
shared_name!Adam/batch_normalization/beta/v
?
3Adam/batch_normalization/beta/v/Read/ReadVariableOpReadVariableOpAdam/batch_normalization/beta/v*
_output_shapes
:*
dtype0
?
Adam/conv3d_1/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*'
shared_nameAdam/conv3d_1/kernel/v
?
*Adam/conv3d_1/kernel/v/Read/ReadVariableOpReadVariableOpAdam/conv3d_1/kernel/v**
_output_shapes
:*
dtype0
?
Adam/conv3d_1/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*%
shared_nameAdam/conv3d_1/bias/v
y
(Adam/conv3d_1/bias/v/Read/ReadVariableOpReadVariableOpAdam/conv3d_1/bias/v*
_output_shapes
:*
dtype0
?
"Adam/batch_normalization_1/gamma/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*3
shared_name$"Adam/batch_normalization_1/gamma/v
?
6Adam/batch_normalization_1/gamma/v/Read/ReadVariableOpReadVariableOp"Adam/batch_normalization_1/gamma/v*
_output_shapes
:*
dtype0
?
!Adam/batch_normalization_1/beta/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*2
shared_name#!Adam/batch_normalization_1/beta/v
?
5Adam/batch_normalization_1/beta/v/Read/ReadVariableOpReadVariableOp!Adam/batch_normalization_1/beta/v*
_output_shapes
:*
dtype0
?
Adam/conv3d_2/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*'
shared_nameAdam/conv3d_2/kernel/v
?
*Adam/conv3d_2/kernel/v/Read/ReadVariableOpReadVariableOpAdam/conv3d_2/kernel/v**
_output_shapes
:*
dtype0
?
Adam/conv3d_2/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*%
shared_nameAdam/conv3d_2/bias/v
y
(Adam/conv3d_2/bias/v/Read/ReadVariableOpReadVariableOpAdam/conv3d_2/bias/v*
_output_shapes
:*
dtype0
?
"Adam/batch_normalization_2/gamma/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*3
shared_name$"Adam/batch_normalization_2/gamma/v
?
6Adam/batch_normalization_2/gamma/v/Read/ReadVariableOpReadVariableOp"Adam/batch_normalization_2/gamma/v*
_output_shapes
:*
dtype0
?
!Adam/batch_normalization_2/beta/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*2
shared_name#!Adam/batch_normalization_2/beta/v
?
5Adam/batch_normalization_2/beta/v/Read/ReadVariableOpReadVariableOp!Adam/batch_normalization_2/beta/v*
_output_shapes
:*
dtype0
?
Adam/conv3d_3/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*'
shared_nameAdam/conv3d_3/kernel/v
?
*Adam/conv3d_3/kernel/v/Read/ReadVariableOpReadVariableOpAdam/conv3d_3/kernel/v**
_output_shapes
:*
dtype0
?
Adam/conv3d_3/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*%
shared_nameAdam/conv3d_3/bias/v
y
(Adam/conv3d_3/bias/v/Read/ReadVariableOpReadVariableOpAdam/conv3d_3/bias/v*
_output_shapes
:*
dtype0
?
"Adam/batch_normalization_3/gamma/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*3
shared_name$"Adam/batch_normalization_3/gamma/v
?
6Adam/batch_normalization_3/gamma/v/Read/ReadVariableOpReadVariableOp"Adam/batch_normalization_3/gamma/v*
_output_shapes
:*
dtype0
?
!Adam/batch_normalization_3/beta/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*2
shared_name#!Adam/batch_normalization_3/beta/v
?
5Adam/batch_normalization_3/beta/v/Read/ReadVariableOpReadVariableOp!Adam/batch_normalization_3/beta/v*
_output_shapes
:*
dtype0
?
Adam/layer1/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
?
?*%
shared_nameAdam/layer1/kernel/v

(Adam/layer1/kernel/v/Read/ReadVariableOpReadVariableOpAdam/layer1/kernel/v* 
_output_shapes
:
?
?*
dtype0
}
Adam/layer1/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:?*#
shared_nameAdam/layer1/bias/v
v
&Adam/layer1/bias/v/Read/ReadVariableOpReadVariableOpAdam/layer1/bias/v*
_output_shapes	
:?*
dtype0
?
Adam/dense/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:	?d*$
shared_nameAdam/dense/kernel/v
|
'Adam/dense/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense/kernel/v*
_output_shapes
:	?d*
dtype0
z
Adam/dense/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:d*"
shared_nameAdam/dense/bias/v
s
%Adam/dense/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense/bias/v*
_output_shapes
:d*
dtype0
?
Adam/dense_1/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:d
*&
shared_nameAdam/dense_1/kernel/v

)Adam/dense_1/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense_1/kernel/v*
_output_shapes

:d
*
dtype0
~
Adam/dense_1/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*$
shared_nameAdam/dense_1/bias/v
w
'Adam/dense_1/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense_1/bias/v*
_output_shapes
:
*
dtype0

NoOpNoOp
??
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*Ř
value??B?? B??
?
layer-0
layer-1
layer_with_weights-0
layer-2
layer-3
layer-4
layer-5
layer_with_weights-1
layer-6
	optimizer
	trainable_variables

regularization_losses
	variables
	keras_api

signatures
 
 
?
layer_with_weights-0
layer-0
layer_with_weights-1
layer-1
layer_with_weights-2
layer-2
trainable_variables
regularization_losses
	variables
	keras_api
R
trainable_variables
regularization_losses
	variables
	keras_api
 
R
trainable_variables
regularization_losses
	variables
	keras_api
h

kernel
bias
trainable_variables
 regularization_losses
!	variables
"	keras_api
?
#iter

$beta_1

%beta_2
	&decay
'learning_ratem?m?(m?)m?*m?+m?,m?-m?.m?/m?0m?1m?2m?3m?4m?5m?6m?7m?8m?9m?:m?;m?<m?=m?v?v?(v?)v?*v?+v?,v?-v?.v?/v?0v?1v?2v?3v?4v?5v?6v?7v?8v?9v?:v?;v?<v?=v?
?
(0
)1
*2
+3
,4
-5
.6
/7
08
19
210
311
412
513
614
715
816
917
:18
;19
<20
=21
22
23
 
?
(0
)1
*2
+3
>4
?5
,6
-7
.8
/9
@10
A11
012
113
214
315
B16
C17
418
519
620
721
D22
E23
824
925
:26
;27
<28
=29
30
31
?
	trainable_variables

regularization_losses
	variables

Flayers
Glayer_metrics
Hmetrics
Inon_trainable_variables
Jlayer_regularization_losses
 
?
Klayer-0
Llayer_with_weights-0
Llayer-1
Mlayer_with_weights-1
Mlayer-2
Nlayer_with_weights-2
Nlayer-3
Olayer_with_weights-3
Olayer-4
Player_with_weights-4
Player-5
Qlayer_with_weights-5
Qlayer-6
Rlayer_with_weights-6
Rlayer-7
Slayer_with_weights-7
Slayer-8
Tlayer-9
Ulayer-10
Vlayer-11
Wlayer_with_weights-8
Wlayer-12
Xtrainable_variables
Yregularization_losses
Z	variables
[	keras_api
h

:kernel
;bias
\trainable_variables
]regularization_losses
^	variables
_	keras_api
h

<kernel
=bias
`trainable_variables
aregularization_losses
b	variables
c	keras_api
?
(0
)1
*2
+3
,4
-5
.6
/7
08
19
210
311
412
513
614
715
816
917
:18
;19
<20
=21
 
?
(0
)1
*2
+3
>4
?5
,6
-7
.8
/9
@10
A11
012
113
214
315
B16
C17
418
519
620
721
D22
E23
824
925
:26
;27
<28
=29
?
trainable_variables
regularization_losses
	variables

dlayers
elayer_metrics
fmetrics
gnon_trainable_variables
hlayer_regularization_losses
 
 
 
?
trainable_variables
regularization_losses
	variables
ilayer_metrics

jlayers
kmetrics
lnon_trainable_variables
mlayer_regularization_losses
 
 
 
?
trainable_variables
regularization_losses
	variables
nlayer_metrics

olayers
pmetrics
qnon_trainable_variables
rlayer_regularization_losses
ZX
VARIABLE_VALUEdense_2/kernel6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUE
VT
VARIABLE_VALUEdense_2/bias4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUE

0
1
 

0
1
?
trainable_variables
 regularization_losses
!	variables
slayer_metrics

tlayers
umetrics
vnon_trainable_variables
wlayer_regularization_losses
HF
VARIABLE_VALUE	Adam/iter)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEAdam/beta_1+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEAdam/beta_2+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUE
JH
VARIABLE_VALUE
Adam/decay*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUE
ZX
VARIABLE_VALUEAdam/learning_rate2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE
SQ
VARIABLE_VALUEconv3d/kernel0trainable_variables/0/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEconv3d/bias0trainable_variables/1/.ATTRIBUTES/VARIABLE_VALUE
_]
VARIABLE_VALUEbatch_normalization/gamma0trainable_variables/2/.ATTRIBUTES/VARIABLE_VALUE
^\
VARIABLE_VALUEbatch_normalization/beta0trainable_variables/3/.ATTRIBUTES/VARIABLE_VALUE
US
VARIABLE_VALUEconv3d_1/kernel0trainable_variables/4/.ATTRIBUTES/VARIABLE_VALUE
SQ
VARIABLE_VALUEconv3d_1/bias0trainable_variables/5/.ATTRIBUTES/VARIABLE_VALUE
a_
VARIABLE_VALUEbatch_normalization_1/gamma0trainable_variables/6/.ATTRIBUTES/VARIABLE_VALUE
`^
VARIABLE_VALUEbatch_normalization_1/beta0trainable_variables/7/.ATTRIBUTES/VARIABLE_VALUE
US
VARIABLE_VALUEconv3d_2/kernel0trainable_variables/8/.ATTRIBUTES/VARIABLE_VALUE
SQ
VARIABLE_VALUEconv3d_2/bias0trainable_variables/9/.ATTRIBUTES/VARIABLE_VALUE
b`
VARIABLE_VALUEbatch_normalization_2/gamma1trainable_variables/10/.ATTRIBUTES/VARIABLE_VALUE
a_
VARIABLE_VALUEbatch_normalization_2/beta1trainable_variables/11/.ATTRIBUTES/VARIABLE_VALUE
VT
VARIABLE_VALUEconv3d_3/kernel1trainable_variables/12/.ATTRIBUTES/VARIABLE_VALUE
TR
VARIABLE_VALUEconv3d_3/bias1trainable_variables/13/.ATTRIBUTES/VARIABLE_VALUE
b`
VARIABLE_VALUEbatch_normalization_3/gamma1trainable_variables/14/.ATTRIBUTES/VARIABLE_VALUE
a_
VARIABLE_VALUEbatch_normalization_3/beta1trainable_variables/15/.ATTRIBUTES/VARIABLE_VALUE
TR
VARIABLE_VALUElayer1/kernel1trainable_variables/16/.ATTRIBUTES/VARIABLE_VALUE
RP
VARIABLE_VALUElayer1/bias1trainable_variables/17/.ATTRIBUTES/VARIABLE_VALUE
SQ
VARIABLE_VALUEdense/kernel1trainable_variables/18/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUE
dense/bias1trainable_variables/19/.ATTRIBUTES/VARIABLE_VALUE
US
VARIABLE_VALUEdense_1/kernel1trainable_variables/20/.ATTRIBUTES/VARIABLE_VALUE
SQ
VARIABLE_VALUEdense_1/bias1trainable_variables/21/.ATTRIBUTES/VARIABLE_VALUE
[Y
VARIABLE_VALUEbatch_normalization/moving_mean&variables/4/.ATTRIBUTES/VARIABLE_VALUE
_]
VARIABLE_VALUE#batch_normalization/moving_variance&variables/5/.ATTRIBUTES/VARIABLE_VALUE
^\
VARIABLE_VALUE!batch_normalization_1/moving_mean'variables/10/.ATTRIBUTES/VARIABLE_VALUE
b`
VARIABLE_VALUE%batch_normalization_1/moving_variance'variables/11/.ATTRIBUTES/VARIABLE_VALUE
^\
VARIABLE_VALUE!batch_normalization_2/moving_mean'variables/16/.ATTRIBUTES/VARIABLE_VALUE
b`
VARIABLE_VALUE%batch_normalization_2/moving_variance'variables/17/.ATTRIBUTES/VARIABLE_VALUE
^\
VARIABLE_VALUE!batch_normalization_3/moving_mean'variables/22/.ATTRIBUTES/VARIABLE_VALUE
b`
VARIABLE_VALUE%batch_normalization_3/moving_variance'variables/23/.ATTRIBUTES/VARIABLE_VALUE
1
0
1
2
3
4
5
6
 

x0
8
>0
?1
@2
A3
B4
C5
D6
E7
 
%
#y_self_saveable_object_factories
?

(kernel
)bias
#z_self_saveable_object_factories
{trainable_variables
|regularization_losses
}	variables
~	keras_api
?
axis
	*gamma
+beta
>moving_mean
?moving_variance
$?_self_saveable_object_factories
?trainable_variables
?regularization_losses
?	variables
?	keras_api
?

,kernel
-bias
$?_self_saveable_object_factories
?trainable_variables
?regularization_losses
?	variables
?	keras_api
?
	?axis
	.gamma
/beta
@moving_mean
Amoving_variance
$?_self_saveable_object_factories
?trainable_variables
?regularization_losses
?	variables
?	keras_api
?

0kernel
1bias
$?_self_saveable_object_factories
?trainable_variables
?regularization_losses
?	variables
?	keras_api
?
	?axis
	2gamma
3beta
Bmoving_mean
Cmoving_variance
$?_self_saveable_object_factories
?trainable_variables
?regularization_losses
?	variables
?	keras_api
?

4kernel
5bias
$?_self_saveable_object_factories
?trainable_variables
?regularization_losses
?	variables
?	keras_api
?
	?axis
	6gamma
7beta
Dmoving_mean
Emoving_variance
$?_self_saveable_object_factories
?trainable_variables
?regularization_losses
?	variables
?	keras_api
|
$?_self_saveable_object_factories
?trainable_variables
?regularization_losses
?	variables
?	keras_api
|
$?_self_saveable_object_factories
?trainable_variables
?regularization_losses
?	variables
?	keras_api
|
$?_self_saveable_object_factories
?trainable_variables
?regularization_losses
?	variables
?	keras_api
?

8kernel
9bias
$?_self_saveable_object_factories
?trainable_variables
?regularization_losses
?	variables
?	keras_api
?
(0
)1
*2
+3
,4
-5
.6
/7
08
19
210
311
412
513
614
715
816
917
 
?
(0
)1
*2
+3
>4
?5
,6
-7
.8
/9
@10
A11
012
113
214
315
B16
C17
418
519
620
721
D22
E23
824
925
?
Xtrainable_variables
Yregularization_losses
Z	variables
?layers
?layer_metrics
?metrics
?non_trainable_variables
 ?layer_regularization_losses

:0
;1
 

:0
;1
?
\trainable_variables
]regularization_losses
^	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses

<0
=1
 

<0
=1
?
`trainable_variables
aregularization_losses
b	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses

0
1
2
 
 
8
>0
?1
@2
A3
B4
C5
D6
E7
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
8

?total

?count
?	variables
?	keras_api
 
 

(0
)1
 

(0
)1
?
{trainable_variables
|regularization_losses
}	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
 
 

*0
+1
 

*0
+1
>2
?3
?
?trainable_variables
?regularization_losses
?	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
 

,0
-1
 

,0
-1
?
?trainable_variables
?regularization_losses
?	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
 
 

.0
/1
 

.0
/1
@2
A3
?
?trainable_variables
?regularization_losses
?	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
 

00
11
 

00
11
?
?trainable_variables
?regularization_losses
?	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
 
 

20
31
 

20
31
B2
C3
?
?trainable_variables
?regularization_losses
?	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
 

40
51
 

40
51
?
?trainable_variables
?regularization_losses
?	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
 
 

60
71
 

60
71
D2
E3
?
?trainable_variables
?regularization_losses
?	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
 
 
 
 
?
?trainable_variables
?regularization_losses
?	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
 
 
 
 
?
?trainable_variables
?regularization_losses
?	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
 
 
 
 
?
?trainable_variables
?regularization_losses
?	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
 

80
91
 

80
91
?
?trainable_variables
?regularization_losses
?	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
^
K0
L1
M2
N3
O4
P5
Q6
R7
S8
T9
U10
V11
W12
 
 
8
>0
?1
@2
A3
B4
C5
D6
E7
 
 
 
 
 
 
 
 
 
 
 
OM
VARIABLE_VALUEtotal4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE
OM
VARIABLE_VALUEcount4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE

?0
?1

?	variables
 
 
 
 
 
 
 
 

>0
?1
 
 
 
 
 
 
 
 
 

@0
A1
 
 
 
 
 
 
 
 
 

B0
C1
 
 
 
 
 
 
 
 
 

D0
E1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
}{
VARIABLE_VALUEAdam/dense_2/kernel/mRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/dense_2/bias/mPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
vt
VARIABLE_VALUEAdam/conv3d/kernel/mLtrainable_variables/0/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
tr
VARIABLE_VALUEAdam/conv3d/bias/mLtrainable_variables/1/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
??
VARIABLE_VALUE Adam/batch_normalization/gamma/mLtrainable_variables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
?
VARIABLE_VALUEAdam/batch_normalization/beta/mLtrainable_variables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/conv3d_1/kernel/mLtrainable_variables/4/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
vt
VARIABLE_VALUEAdam/conv3d_1/bias/mLtrainable_variables/5/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
??
VARIABLE_VALUE"Adam/batch_normalization_1/gamma/mLtrainable_variables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
??
VARIABLE_VALUE!Adam/batch_normalization_1/beta/mLtrainable_variables/7/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/conv3d_2/kernel/mLtrainable_variables/8/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
vt
VARIABLE_VALUEAdam/conv3d_2/bias/mLtrainable_variables/9/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
??
VARIABLE_VALUE"Adam/batch_normalization_2/gamma/mMtrainable_variables/10/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
??
VARIABLE_VALUE!Adam/batch_normalization_2/beta/mMtrainable_variables/11/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/conv3d_3/kernel/mMtrainable_variables/12/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
wu
VARIABLE_VALUEAdam/conv3d_3/bias/mMtrainable_variables/13/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
??
VARIABLE_VALUE"Adam/batch_normalization_3/gamma/mMtrainable_variables/14/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
??
VARIABLE_VALUE!Adam/batch_normalization_3/beta/mMtrainable_variables/15/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
wu
VARIABLE_VALUEAdam/layer1/kernel/mMtrainable_variables/16/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
us
VARIABLE_VALUEAdam/layer1/bias/mMtrainable_variables/17/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
vt
VARIABLE_VALUEAdam/dense/kernel/mMtrainable_variables/18/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
tr
VARIABLE_VALUEAdam/dense/bias/mMtrainable_variables/19/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/dense_1/kernel/mMtrainable_variables/20/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
vt
VARIABLE_VALUEAdam/dense_1/bias/mMtrainable_variables/21/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
}{
VARIABLE_VALUEAdam/dense_2/kernel/vRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/dense_2/bias/vPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
vt
VARIABLE_VALUEAdam/conv3d/kernel/vLtrainable_variables/0/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
tr
VARIABLE_VALUEAdam/conv3d/bias/vLtrainable_variables/1/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
??
VARIABLE_VALUE Adam/batch_normalization/gamma/vLtrainable_variables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
?
VARIABLE_VALUEAdam/batch_normalization/beta/vLtrainable_variables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/conv3d_1/kernel/vLtrainable_variables/4/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
vt
VARIABLE_VALUEAdam/conv3d_1/bias/vLtrainable_variables/5/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
??
VARIABLE_VALUE"Adam/batch_normalization_1/gamma/vLtrainable_variables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
??
VARIABLE_VALUE!Adam/batch_normalization_1/beta/vLtrainable_variables/7/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/conv3d_2/kernel/vLtrainable_variables/8/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
vt
VARIABLE_VALUEAdam/conv3d_2/bias/vLtrainable_variables/9/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
??
VARIABLE_VALUE"Adam/batch_normalization_2/gamma/vMtrainable_variables/10/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
??
VARIABLE_VALUE!Adam/batch_normalization_2/beta/vMtrainable_variables/11/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/conv3d_3/kernel/vMtrainable_variables/12/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
wu
VARIABLE_VALUEAdam/conv3d_3/bias/vMtrainable_variables/13/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
??
VARIABLE_VALUE"Adam/batch_normalization_3/gamma/vMtrainable_variables/14/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
??
VARIABLE_VALUE!Adam/batch_normalization_3/beta/vMtrainable_variables/15/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
wu
VARIABLE_VALUEAdam/layer1/kernel/vMtrainable_variables/16/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
us
VARIABLE_VALUEAdam/layer1/bias/vMtrainable_variables/17/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
vt
VARIABLE_VALUEAdam/dense/kernel/vMtrainable_variables/18/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
tr
VARIABLE_VALUEAdam/dense/bias/vMtrainable_variables/19/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
xv
VARIABLE_VALUEAdam/dense_1/kernel/vMtrainable_variables/20/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
vt
VARIABLE_VALUEAdam/dense_1/bias/vMtrainable_variables/21/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
?
serving_default_input_1Placeholder*4
_output_shapes"
 :??????????*
dtype0*)
shape :??????????
?
serving_default_input_2Placeholder*4
_output_shapes"
 :??????????*
dtype0*)
shape :??????????
z
serving_default_input_3Placeholder*'
_output_shapes
:?????????
*
dtype0*
shape:?????????

?	
StatefulPartitionedCallStatefulPartitionedCallserving_default_input_1serving_default_input_2serving_default_input_3conv3d/kernelconv3d/bias#batch_normalization/moving_variancebatch_normalization/gammabatch_normalization/moving_meanbatch_normalization/betaconv3d_1/kernelconv3d_1/bias%batch_normalization_1/moving_variancebatch_normalization_1/gamma!batch_normalization_1/moving_meanbatch_normalization_1/betaconv3d_2/kernelconv3d_2/bias%batch_normalization_2/moving_variancebatch_normalization_2/gamma!batch_normalization_2/moving_meanbatch_normalization_2/betaconv3d_3/kernelconv3d_3/bias%batch_normalization_3/moving_variancebatch_normalization_3/gamma!batch_normalization_3/moving_meanbatch_normalization_3/betalayer1/kernellayer1/biasdense/kernel
dense/biasdense_1/kerneldense_1/biasdense_2/kerneldense_2/bias*.
Tin'
%2#*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????*B
_read_only_resource_inputs$
" 	
 !"*0
config_proto 

CPU

GPU2*0J 8? *.
f)R'
%__inference_signature_wrapper_1225148
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
? 
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filename"dense_2/kernel/Read/ReadVariableOp dense_2/bias/Read/ReadVariableOpAdam/iter/Read/ReadVariableOpAdam/beta_1/Read/ReadVariableOpAdam/beta_2/Read/ReadVariableOpAdam/decay/Read/ReadVariableOp&Adam/learning_rate/Read/ReadVariableOp!conv3d/kernel/Read/ReadVariableOpconv3d/bias/Read/ReadVariableOp-batch_normalization/gamma/Read/ReadVariableOp,batch_normalization/beta/Read/ReadVariableOp#conv3d_1/kernel/Read/ReadVariableOp!conv3d_1/bias/Read/ReadVariableOp/batch_normalization_1/gamma/Read/ReadVariableOp.batch_normalization_1/beta/Read/ReadVariableOp#conv3d_2/kernel/Read/ReadVariableOp!conv3d_2/bias/Read/ReadVariableOp/batch_normalization_2/gamma/Read/ReadVariableOp.batch_normalization_2/beta/Read/ReadVariableOp#conv3d_3/kernel/Read/ReadVariableOp!conv3d_3/bias/Read/ReadVariableOp/batch_normalization_3/gamma/Read/ReadVariableOp.batch_normalization_3/beta/Read/ReadVariableOp!layer1/kernel/Read/ReadVariableOplayer1/bias/Read/ReadVariableOp dense/kernel/Read/ReadVariableOpdense/bias/Read/ReadVariableOp"dense_1/kernel/Read/ReadVariableOp dense_1/bias/Read/ReadVariableOp3batch_normalization/moving_mean/Read/ReadVariableOp7batch_normalization/moving_variance/Read/ReadVariableOp5batch_normalization_1/moving_mean/Read/ReadVariableOp9batch_normalization_1/moving_variance/Read/ReadVariableOp5batch_normalization_2/moving_mean/Read/ReadVariableOp9batch_normalization_2/moving_variance/Read/ReadVariableOp5batch_normalization_3/moving_mean/Read/ReadVariableOp9batch_normalization_3/moving_variance/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOp)Adam/dense_2/kernel/m/Read/ReadVariableOp'Adam/dense_2/bias/m/Read/ReadVariableOp(Adam/conv3d/kernel/m/Read/ReadVariableOp&Adam/conv3d/bias/m/Read/ReadVariableOp4Adam/batch_normalization/gamma/m/Read/ReadVariableOp3Adam/batch_normalization/beta/m/Read/ReadVariableOp*Adam/conv3d_1/kernel/m/Read/ReadVariableOp(Adam/conv3d_1/bias/m/Read/ReadVariableOp6Adam/batch_normalization_1/gamma/m/Read/ReadVariableOp5Adam/batch_normalization_1/beta/m/Read/ReadVariableOp*Adam/conv3d_2/kernel/m/Read/ReadVariableOp(Adam/conv3d_2/bias/m/Read/ReadVariableOp6Adam/batch_normalization_2/gamma/m/Read/ReadVariableOp5Adam/batch_normalization_2/beta/m/Read/ReadVariableOp*Adam/conv3d_3/kernel/m/Read/ReadVariableOp(Adam/conv3d_3/bias/m/Read/ReadVariableOp6Adam/batch_normalization_3/gamma/m/Read/ReadVariableOp5Adam/batch_normalization_3/beta/m/Read/ReadVariableOp(Adam/layer1/kernel/m/Read/ReadVariableOp&Adam/layer1/bias/m/Read/ReadVariableOp'Adam/dense/kernel/m/Read/ReadVariableOp%Adam/dense/bias/m/Read/ReadVariableOp)Adam/dense_1/kernel/m/Read/ReadVariableOp'Adam/dense_1/bias/m/Read/ReadVariableOp)Adam/dense_2/kernel/v/Read/ReadVariableOp'Adam/dense_2/bias/v/Read/ReadVariableOp(Adam/conv3d/kernel/v/Read/ReadVariableOp&Adam/conv3d/bias/v/Read/ReadVariableOp4Adam/batch_normalization/gamma/v/Read/ReadVariableOp3Adam/batch_normalization/beta/v/Read/ReadVariableOp*Adam/conv3d_1/kernel/v/Read/ReadVariableOp(Adam/conv3d_1/bias/v/Read/ReadVariableOp6Adam/batch_normalization_1/gamma/v/Read/ReadVariableOp5Adam/batch_normalization_1/beta/v/Read/ReadVariableOp*Adam/conv3d_2/kernel/v/Read/ReadVariableOp(Adam/conv3d_2/bias/v/Read/ReadVariableOp6Adam/batch_normalization_2/gamma/v/Read/ReadVariableOp5Adam/batch_normalization_2/beta/v/Read/ReadVariableOp*Adam/conv3d_3/kernel/v/Read/ReadVariableOp(Adam/conv3d_3/bias/v/Read/ReadVariableOp6Adam/batch_normalization_3/gamma/v/Read/ReadVariableOp5Adam/batch_normalization_3/beta/v/Read/ReadVariableOp(Adam/layer1/kernel/v/Read/ReadVariableOp&Adam/layer1/bias/v/Read/ReadVariableOp'Adam/dense/kernel/v/Read/ReadVariableOp%Adam/dense/bias/v/Read/ReadVariableOp)Adam/dense_1/kernel/v/Read/ReadVariableOp'Adam/dense_1/bias/v/Read/ReadVariableOpConst*d
Tin]
[2Y	*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *)
f$R"
 __inference__traced_save_1227989
?
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenamedense_2/kerneldense_2/bias	Adam/iterAdam/beta_1Adam/beta_2
Adam/decayAdam/learning_rateconv3d/kernelconv3d/biasbatch_normalization/gammabatch_normalization/betaconv3d_1/kernelconv3d_1/biasbatch_normalization_1/gammabatch_normalization_1/betaconv3d_2/kernelconv3d_2/biasbatch_normalization_2/gammabatch_normalization_2/betaconv3d_3/kernelconv3d_3/biasbatch_normalization_3/gammabatch_normalization_3/betalayer1/kernellayer1/biasdense/kernel
dense/biasdense_1/kerneldense_1/biasbatch_normalization/moving_mean#batch_normalization/moving_variance!batch_normalization_1/moving_mean%batch_normalization_1/moving_variance!batch_normalization_2/moving_mean%batch_normalization_2/moving_variance!batch_normalization_3/moving_mean%batch_normalization_3/moving_variancetotalcountAdam/dense_2/kernel/mAdam/dense_2/bias/mAdam/conv3d/kernel/mAdam/conv3d/bias/m Adam/batch_normalization/gamma/mAdam/batch_normalization/beta/mAdam/conv3d_1/kernel/mAdam/conv3d_1/bias/m"Adam/batch_normalization_1/gamma/m!Adam/batch_normalization_1/beta/mAdam/conv3d_2/kernel/mAdam/conv3d_2/bias/m"Adam/batch_normalization_2/gamma/m!Adam/batch_normalization_2/beta/mAdam/conv3d_3/kernel/mAdam/conv3d_3/bias/m"Adam/batch_normalization_3/gamma/m!Adam/batch_normalization_3/beta/mAdam/layer1/kernel/mAdam/layer1/bias/mAdam/dense/kernel/mAdam/dense/bias/mAdam/dense_1/kernel/mAdam/dense_1/bias/mAdam/dense_2/kernel/vAdam/dense_2/bias/vAdam/conv3d/kernel/vAdam/conv3d/bias/v Adam/batch_normalization/gamma/vAdam/batch_normalization/beta/vAdam/conv3d_1/kernel/vAdam/conv3d_1/bias/v"Adam/batch_normalization_1/gamma/v!Adam/batch_normalization_1/beta/vAdam/conv3d_2/kernel/vAdam/conv3d_2/bias/v"Adam/batch_normalization_2/gamma/v!Adam/batch_normalization_2/beta/vAdam/conv3d_3/kernel/vAdam/conv3d_3/bias/v"Adam/batch_normalization_3/gamma/v!Adam/batch_normalization_3/beta/vAdam/layer1/kernel/vAdam/layer1/bias/vAdam/dense/kernel/vAdam/dense/bias/vAdam/dense_1/kernel/vAdam/dense_1/bias/v*c
Tin\
Z2X*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *,
f'R%
#__inference__traced_restore_1228260??+
?
?
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_1227353

inputs%
!batchnorm_readvariableop_resource)
%batchnorm_mul_readvariableop_resource'
#batchnorm_readvariableop_1_resource'
#batchnorm_readvariableop_2_resource
identity??
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/mul_1?
batchnorm/ReadVariableOp_1ReadVariableOp#batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_1?
batchnorm/mul_2Mul"batchnorm/ReadVariableOp_1:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOp_2ReadVariableOp#batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_2?
batchnorm/subSub"batchnorm/ReadVariableOp_2:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????:::::v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
?,
?
G__inference_sequential_layer_call_and_return_conditional_losses_1224186

inputs
functional_1_1224110
functional_1_1224112
functional_1_1224114
functional_1_1224116
functional_1_1224118
functional_1_1224120
functional_1_1224122
functional_1_1224124
functional_1_1224126
functional_1_1224128
functional_1_1224130
functional_1_1224132
functional_1_1224134
functional_1_1224136
functional_1_1224138
functional_1_1224140
functional_1_1224142
functional_1_1224144
functional_1_1224146
functional_1_1224148
functional_1_1224150
functional_1_1224152
functional_1_1224154
functional_1_1224156
functional_1_1224158
functional_1_1224160
dense_1224163
dense_1224165
dense_1_1224168
dense_1_1224170
identity??dense/StatefulPartitionedCall?dense_1/StatefulPartitionedCall?$functional_1/StatefulPartitionedCall?
$functional_1/StatefulPartitionedCallStatefulPartitionedCallinputsfunctional_1_1224110functional_1_1224112functional_1_1224114functional_1_1224116functional_1_1224118functional_1_1224120functional_1_1224122functional_1_1224124functional_1_1224126functional_1_1224128functional_1_1224130functional_1_1224132functional_1_1224134functional_1_1224136functional_1_1224138functional_1_1224140functional_1_1224142functional_1_1224144functional_1_1224146functional_1_1224148functional_1_1224150functional_1_1224152functional_1_1224154functional_1_1224156functional_1_1224158functional_1_1224160*&
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*<
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8? *R
fMRK
I__inference_functional_1_layer_call_and_return_conditional_losses_12235762&
$functional_1/StatefulPartitionedCall?
dense/StatefulPartitionedCallStatefulPartitionedCall-functional_1/StatefulPartitionedCall:output:0dense_1224163dense_1224165*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????d*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *K
fFRD
B__inference_dense_layer_call_and_return_conditional_losses_12238192
dense/StatefulPartitionedCall?
dense_1/StatefulPartitionedCallStatefulPartitionedCall&dense/StatefulPartitionedCall:output:0dense_1_1224168dense_1_1224170*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_1_layer_call_and_return_conditional_losses_12238522!
dense_1/StatefulPartitionedCall?
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_1224163*
_output_shapes
:	?d*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOp?
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?d2!
dense/kernel/Regularizer/Square?
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const?
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/Sum?
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2 
dense/kernel/Regularizer/mul/x?
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/mul?
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_1224168*
_output_shapes

:d
*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp?
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d
2#
!dense_1/kernel/Regularizer/Square?
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const?
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/Sum?
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2"
 dense_1/kernel/Regularizer/mul/x?
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/mul?
IdentityIdentity(dense_1/StatefulPartitionedCall:output:0^dense/StatefulPartitionedCall ^dense_1/StatefulPartitionedCall%^functional_1/StatefulPartitionedCall*
T0*'
_output_shapes
:?????????
2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????::::::::::::::::::::::::::::::2>
dense/StatefulPartitionedCalldense/StatefulPartitionedCall2B
dense_1/StatefulPartitionedCalldense_1/StatefulPartitionedCall2L
$functional_1/StatefulPartitionedCall$functional_1/StatefulPartitionedCall:\ X
4
_output_shapes"
 :??????????
 
_user_specified_nameinputs
?
?
7__inference_batch_normalization_2_layer_call_fn_1227461

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_12230852
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????::::22
StatefulPartitionedCallStatefulPartitionedCall:[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
??
?
I__inference_functional_1_layer_call_and_return_conditional_losses_1226710

inputs)
%conv3d_conv3d_readvariableop_resource*
&conv3d_biasadd_readvariableop_resource9
5batch_normalization_batchnorm_readvariableop_resource=
9batch_normalization_batchnorm_mul_readvariableop_resource;
7batch_normalization_batchnorm_readvariableop_1_resource;
7batch_normalization_batchnorm_readvariableop_2_resource+
'conv3d_1_conv3d_readvariableop_resource,
(conv3d_1_biasadd_readvariableop_resource;
7batch_normalization_1_batchnorm_readvariableop_resource?
;batch_normalization_1_batchnorm_mul_readvariableop_resource=
9batch_normalization_1_batchnorm_readvariableop_1_resource=
9batch_normalization_1_batchnorm_readvariableop_2_resource+
'conv3d_2_conv3d_readvariableop_resource,
(conv3d_2_biasadd_readvariableop_resource;
7batch_normalization_2_batchnorm_readvariableop_resource?
;batch_normalization_2_batchnorm_mul_readvariableop_resource=
9batch_normalization_2_batchnorm_readvariableop_1_resource=
9batch_normalization_2_batchnorm_readvariableop_2_resource+
'conv3d_3_conv3d_readvariableop_resource,
(conv3d_3_biasadd_readvariableop_resource;
7batch_normalization_3_batchnorm_readvariableop_resource?
;batch_normalization_3_batchnorm_mul_readvariableop_resource=
9batch_normalization_3_batchnorm_readvariableop_1_resource=
9batch_normalization_3_batchnorm_readvariableop_2_resource)
%layer1_matmul_readvariableop_resource*
&layer1_biasadd_readvariableop_resource
identity??
conv3d/Conv3D/ReadVariableOpReadVariableOp%conv3d_conv3d_readvariableop_resource*+
_output_shapes
:?*
dtype02
conv3d/Conv3D/ReadVariableOp?
conv3d/Conv3DConv3Dinputs$conv3d/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
conv3d/Conv3D?
conv3d/BiasAdd/ReadVariableOpReadVariableOp&conv3d_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02
conv3d/BiasAdd/ReadVariableOp?
conv3d/BiasAddBiasAddconv3d/Conv3D:output:0%conv3d/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2
conv3d/BiasAdd?
,batch_normalization/batchnorm/ReadVariableOpReadVariableOp5batch_normalization_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02.
,batch_normalization/batchnorm/ReadVariableOp?
#batch_normalization/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2%
#batch_normalization/batchnorm/add/y?
!batch_normalization/batchnorm/addAddV24batch_normalization/batchnorm/ReadVariableOp:value:0,batch_normalization/batchnorm/add/y:output:0*
T0*
_output_shapes
:2#
!batch_normalization/batchnorm/add?
#batch_normalization/batchnorm/RsqrtRsqrt%batch_normalization/batchnorm/add:z:0*
T0*
_output_shapes
:2%
#batch_normalization/batchnorm/Rsqrt?
0batch_normalization/batchnorm/mul/ReadVariableOpReadVariableOp9batch_normalization_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype022
0batch_normalization/batchnorm/mul/ReadVariableOp?
!batch_normalization/batchnorm/mulMul'batch_normalization/batchnorm/Rsqrt:y:08batch_normalization/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2#
!batch_normalization/batchnorm/mul?
#batch_normalization/batchnorm/mul_1Mulconv3d/BiasAdd:output:0%batch_normalization/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2%
#batch_normalization/batchnorm/mul_1?
.batch_normalization/batchnorm/ReadVariableOp_1ReadVariableOp7batch_normalization_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype020
.batch_normalization/batchnorm/ReadVariableOp_1?
#batch_normalization/batchnorm/mul_2Mul6batch_normalization/batchnorm/ReadVariableOp_1:value:0%batch_normalization/batchnorm/mul:z:0*
T0*
_output_shapes
:2%
#batch_normalization/batchnorm/mul_2?
.batch_normalization/batchnorm/ReadVariableOp_2ReadVariableOp7batch_normalization_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype020
.batch_normalization/batchnorm/ReadVariableOp_2?
!batch_normalization/batchnorm/subSub6batch_normalization/batchnorm/ReadVariableOp_2:value:0'batch_normalization/batchnorm/mul_2:z:0*
T0*
_output_shapes
:2#
!batch_normalization/batchnorm/sub?
#batch_normalization/batchnorm/add_1AddV2'batch_normalization/batchnorm/mul_1:z:0%batch_normalization/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2%
#batch_normalization/batchnorm/add_1?
conv3d_1/Conv3D/ReadVariableOpReadVariableOp'conv3d_1_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02 
conv3d_1/Conv3D/ReadVariableOp?
conv3d_1/Conv3DConv3D'batch_normalization/batchnorm/add_1:z:0&conv3d_1/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
conv3d_1/Conv3D?
conv3d_1/BiasAdd/ReadVariableOpReadVariableOp(conv3d_1_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02!
conv3d_1/BiasAdd/ReadVariableOp?
conv3d_1/BiasAddBiasAddconv3d_1/Conv3D:output:0'conv3d_1/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2
conv3d_1/BiasAdd|
conv3d_1/EluEluconv3d_1/BiasAdd:output:0*
T0*3
_output_shapes!
:?????????2
conv3d_1/Elu?
.batch_normalization_1/batchnorm/ReadVariableOpReadVariableOp7batch_normalization_1_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype020
.batch_normalization_1/batchnorm/ReadVariableOp?
%batch_normalization_1/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2'
%batch_normalization_1/batchnorm/add/y?
#batch_normalization_1/batchnorm/addAddV26batch_normalization_1/batchnorm/ReadVariableOp:value:0.batch_normalization_1/batchnorm/add/y:output:0*
T0*
_output_shapes
:2%
#batch_normalization_1/batchnorm/add?
%batch_normalization_1/batchnorm/RsqrtRsqrt'batch_normalization_1/batchnorm/add:z:0*
T0*
_output_shapes
:2'
%batch_normalization_1/batchnorm/Rsqrt?
2batch_normalization_1/batchnorm/mul/ReadVariableOpReadVariableOp;batch_normalization_1_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype024
2batch_normalization_1/batchnorm/mul/ReadVariableOp?
#batch_normalization_1/batchnorm/mulMul)batch_normalization_1/batchnorm/Rsqrt:y:0:batch_normalization_1/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2%
#batch_normalization_1/batchnorm/mul?
%batch_normalization_1/batchnorm/mul_1Mulconv3d_1/Elu:activations:0'batch_normalization_1/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2'
%batch_normalization_1/batchnorm/mul_1?
0batch_normalization_1/batchnorm/ReadVariableOp_1ReadVariableOp9batch_normalization_1_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype022
0batch_normalization_1/batchnorm/ReadVariableOp_1?
%batch_normalization_1/batchnorm/mul_2Mul8batch_normalization_1/batchnorm/ReadVariableOp_1:value:0'batch_normalization_1/batchnorm/mul:z:0*
T0*
_output_shapes
:2'
%batch_normalization_1/batchnorm/mul_2?
0batch_normalization_1/batchnorm/ReadVariableOp_2ReadVariableOp9batch_normalization_1_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype022
0batch_normalization_1/batchnorm/ReadVariableOp_2?
#batch_normalization_1/batchnorm/subSub8batch_normalization_1/batchnorm/ReadVariableOp_2:value:0)batch_normalization_1/batchnorm/mul_2:z:0*
T0*
_output_shapes
:2%
#batch_normalization_1/batchnorm/sub?
%batch_normalization_1/batchnorm/add_1AddV2)batch_normalization_1/batchnorm/mul_1:z:0'batch_normalization_1/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2'
%batch_normalization_1/batchnorm/add_1?
conv3d_2/Conv3D/ReadVariableOpReadVariableOp'conv3d_2_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02 
conv3d_2/Conv3D/ReadVariableOp?
conv3d_2/Conv3DConv3D)batch_normalization_1/batchnorm/add_1:z:0&conv3d_2/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
conv3d_2/Conv3D?
conv3d_2/BiasAdd/ReadVariableOpReadVariableOp(conv3d_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02!
conv3d_2/BiasAdd/ReadVariableOp?
conv3d_2/BiasAddBiasAddconv3d_2/Conv3D:output:0'conv3d_2/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2
conv3d_2/BiasAdd|
conv3d_2/EluEluconv3d_2/BiasAdd:output:0*
T0*3
_output_shapes!
:?????????2
conv3d_2/Elu?
.batch_normalization_2/batchnorm/ReadVariableOpReadVariableOp7batch_normalization_2_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype020
.batch_normalization_2/batchnorm/ReadVariableOp?
%batch_normalization_2/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2'
%batch_normalization_2/batchnorm/add/y?
#batch_normalization_2/batchnorm/addAddV26batch_normalization_2/batchnorm/ReadVariableOp:value:0.batch_normalization_2/batchnorm/add/y:output:0*
T0*
_output_shapes
:2%
#batch_normalization_2/batchnorm/add?
%batch_normalization_2/batchnorm/RsqrtRsqrt'batch_normalization_2/batchnorm/add:z:0*
T0*
_output_shapes
:2'
%batch_normalization_2/batchnorm/Rsqrt?
2batch_normalization_2/batchnorm/mul/ReadVariableOpReadVariableOp;batch_normalization_2_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype024
2batch_normalization_2/batchnorm/mul/ReadVariableOp?
#batch_normalization_2/batchnorm/mulMul)batch_normalization_2/batchnorm/Rsqrt:y:0:batch_normalization_2/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2%
#batch_normalization_2/batchnorm/mul?
%batch_normalization_2/batchnorm/mul_1Mulconv3d_2/Elu:activations:0'batch_normalization_2/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2'
%batch_normalization_2/batchnorm/mul_1?
0batch_normalization_2/batchnorm/ReadVariableOp_1ReadVariableOp9batch_normalization_2_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype022
0batch_normalization_2/batchnorm/ReadVariableOp_1?
%batch_normalization_2/batchnorm/mul_2Mul8batch_normalization_2/batchnorm/ReadVariableOp_1:value:0'batch_normalization_2/batchnorm/mul:z:0*
T0*
_output_shapes
:2'
%batch_normalization_2/batchnorm/mul_2?
0batch_normalization_2/batchnorm/ReadVariableOp_2ReadVariableOp9batch_normalization_2_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype022
0batch_normalization_2/batchnorm/ReadVariableOp_2?
#batch_normalization_2/batchnorm/subSub8batch_normalization_2/batchnorm/ReadVariableOp_2:value:0)batch_normalization_2/batchnorm/mul_2:z:0*
T0*
_output_shapes
:2%
#batch_normalization_2/batchnorm/sub?
%batch_normalization_2/batchnorm/add_1AddV2)batch_normalization_2/batchnorm/mul_1:z:0'batch_normalization_2/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2'
%batch_normalization_2/batchnorm/add_1?
conv3d_3/Conv3D/ReadVariableOpReadVariableOp'conv3d_3_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02 
conv3d_3/Conv3D/ReadVariableOp?
conv3d_3/Conv3DConv3D)batch_normalization_2/batchnorm/add_1:z:0&conv3d_3/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
conv3d_3/Conv3D?
conv3d_3/BiasAdd/ReadVariableOpReadVariableOp(conv3d_3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02!
conv3d_3/BiasAdd/ReadVariableOp?
conv3d_3/BiasAddBiasAddconv3d_3/Conv3D:output:0'conv3d_3/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2
conv3d_3/BiasAdd|
conv3d_3/EluEluconv3d_3/BiasAdd:output:0*
T0*3
_output_shapes!
:?????????2
conv3d_3/Elu?
.batch_normalization_3/batchnorm/ReadVariableOpReadVariableOp7batch_normalization_3_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype020
.batch_normalization_3/batchnorm/ReadVariableOp?
%batch_normalization_3/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2'
%batch_normalization_3/batchnorm/add/y?
#batch_normalization_3/batchnorm/addAddV26batch_normalization_3/batchnorm/ReadVariableOp:value:0.batch_normalization_3/batchnorm/add/y:output:0*
T0*
_output_shapes
:2%
#batch_normalization_3/batchnorm/add?
%batch_normalization_3/batchnorm/RsqrtRsqrt'batch_normalization_3/batchnorm/add:z:0*
T0*
_output_shapes
:2'
%batch_normalization_3/batchnorm/Rsqrt?
2batch_normalization_3/batchnorm/mul/ReadVariableOpReadVariableOp;batch_normalization_3_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype024
2batch_normalization_3/batchnorm/mul/ReadVariableOp?
#batch_normalization_3/batchnorm/mulMul)batch_normalization_3/batchnorm/Rsqrt:y:0:batch_normalization_3/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2%
#batch_normalization_3/batchnorm/mul?
%batch_normalization_3/batchnorm/mul_1Mulconv3d_3/Elu:activations:0'batch_normalization_3/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2'
%batch_normalization_3/batchnorm/mul_1?
0batch_normalization_3/batchnorm/ReadVariableOp_1ReadVariableOp9batch_normalization_3_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype022
0batch_normalization_3/batchnorm/ReadVariableOp_1?
%batch_normalization_3/batchnorm/mul_2Mul8batch_normalization_3/batchnorm/ReadVariableOp_1:value:0'batch_normalization_3/batchnorm/mul:z:0*
T0*
_output_shapes
:2'
%batch_normalization_3/batchnorm/mul_2?
0batch_normalization_3/batchnorm/ReadVariableOp_2ReadVariableOp9batch_normalization_3_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype022
0batch_normalization_3/batchnorm/ReadVariableOp_2?
#batch_normalization_3/batchnorm/subSub8batch_normalization_3/batchnorm/ReadVariableOp_2:value:0)batch_normalization_3/batchnorm/mul_2:z:0*
T0*
_output_shapes
:2%
#batch_normalization_3/batchnorm/sub?
%batch_normalization_3/batchnorm/add_1AddV2)batch_normalization_3/batchnorm/mul_1:z:0'batch_normalization_3/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2'
%batch_normalization_3/batchnorm/add_1?
average_pooling3d/AvgPool3D	AvgPool3D)batch_normalization_3/batchnorm/add_1:z:0*
T0*3
_output_shapes!
:?????????*
ksize	
*
paddingVALID*
strides	
2
average_pooling3d/AvgPool3Do
flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"????   2
flatten/Const?
flatten/ReshapeReshape$average_pooling3d/AvgPool3D:output:0flatten/Const:output:0*
T0*(
_output_shapes
:??????????
2
flatten/Reshape}
dropout/IdentityIdentityflatten/Reshape:output:0*
T0*(
_output_shapes
:??????????
2
dropout/Identity?
layer1/MatMul/ReadVariableOpReadVariableOp%layer1_matmul_readvariableop_resource* 
_output_shapes
:
?
?*
dtype02
layer1/MatMul/ReadVariableOp?
layer1/MatMulMatMuldropout/Identity:output:0$layer1/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2
layer1/MatMul?
layer1/BiasAdd/ReadVariableOpReadVariableOp&layer1_biasadd_readvariableop_resource*
_output_shapes	
:?*
dtype02
layer1/BiasAdd/ReadVariableOp?
layer1/BiasAddBiasAddlayer1/MatMul:product:0%layer1/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2
layer1/BiasAddk

layer1/EluElulayer1/BiasAdd:output:0*
T0*(
_output_shapes
:??????????2

layer1/Elum
IdentityIdentitylayer1/Elu:activations:0*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????:::::::::::::::::::::::::::\ X
4
_output_shapes"
 :??????????
 
_user_specified_nameinputs
?+
?
P__inference_batch_normalization_layer_call_and_return_conditional_losses_1222288

inputs
assignmovingavg_1222263
assignmovingavg_1_1222269)
%batchnorm_mul_readvariableop_resource%
!batchnorm_readvariableop_resource
identity??#AssignMovingAvg/AssignSubVariableOp?%AssignMovingAvg_1/AssignSubVariableOp?
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2 
moments/mean/reduction_indices?
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/mean?
moments/StopGradientStopGradientmoments/mean:output:0*
T0**
_output_shapes
:2
moments/StopGradient?
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
moments/SquaredDifference?
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2$
"moments/variance/reduction_indices?
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/variance?
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze?
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze_1?
AssignMovingAvg/decayConst**
_class 
loc:@AssignMovingAvg/1222263*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg/decay?
AssignMovingAvg/ReadVariableOpReadVariableOpassignmovingavg_1222263*
_output_shapes
:*
dtype02 
AssignMovingAvg/ReadVariableOp?
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0**
_class 
loc:@AssignMovingAvg/1222263*
_output_shapes
:2
AssignMovingAvg/sub?
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0**
_class 
loc:@AssignMovingAvg/1222263*
_output_shapes
:2
AssignMovingAvg/mul?
#AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1222263AssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp**
_class 
loc:@AssignMovingAvg/1222263*
_output_shapes
 *
dtype02%
#AssignMovingAvg/AssignSubVariableOp?
AssignMovingAvg_1/decayConst*,
_class"
 loc:@AssignMovingAvg_1/1222269*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg_1/decay?
 AssignMovingAvg_1/ReadVariableOpReadVariableOpassignmovingavg_1_1222269*
_output_shapes
:*
dtype02"
 AssignMovingAvg_1/ReadVariableOp?
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1222269*
_output_shapes
:2
AssignMovingAvg_1/sub?
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1222269*
_output_shapes
:2
AssignMovingAvg_1/mul?
%AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1_1222269AssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*,
_class"
 loc:@AssignMovingAvg_1/1222269*
_output_shapes
 *
dtype02'
%AssignMovingAvg_1/AssignSubVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/mul_1{
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp?
batchnorm/subSub batchnorm/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0$^AssignMovingAvg/AssignSubVariableOp&^AssignMovingAvg_1/AssignSubVariableOp*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????::::2J
#AssignMovingAvg/AssignSubVariableOp#AssignMovingAvg/AssignSubVariableOp2N
%AssignMovingAvg_1/AssignSubVariableOp%AssignMovingAvg_1/AssignSubVariableOp:v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_1222461

inputs%
!batchnorm_readvariableop_resource)
%batchnorm_mul_readvariableop_resource'
#batchnorm_readvariableop_1_resource'
#batchnorm_readvariableop_2_resource
identity??
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/mul_1?
batchnorm/ReadVariableOp_1ReadVariableOp#batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_1?
batchnorm/mul_2Mul"batchnorm/ReadVariableOp_1:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOp_2ReadVariableOp#batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_2?
batchnorm/subSub"batchnorm/ReadVariableOp_2:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????:::::v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
?
?
B__inference_dense_layer_call_and_return_conditional_losses_1226847

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity??
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	?d*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????d2
MatMul?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:d*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????d2	
BiasAddU
EluEluBiasAdd:output:0*
T0*'
_output_shapes
:?????????d2
Elu?
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	?d*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOp?
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?d2!
dense/kernel/Regularizer/Square?
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const?
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/Sum?
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2 
dense/kernel/Regularizer/mul/x?
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/mule
IdentityIdentityElu:activations:0*
T0*'
_output_shapes
:?????????d2

Identity"
identityIdentity:output:0*/
_input_shapes
:??????????:::P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs
?
?
7__inference_batch_normalization_2_layer_call_fn_1227366

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *N
_output_shapes<
::8????????????????????????????????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_12225682
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????::::22
StatefulPartitionedCallStatefulPartitionedCall:v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
?+
?
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_1222428

inputs
assignmovingavg_1222403
assignmovingavg_1_1222409)
%batchnorm_mul_readvariableop_resource%
!batchnorm_readvariableop_resource
identity??#AssignMovingAvg/AssignSubVariableOp?%AssignMovingAvg_1/AssignSubVariableOp?
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2 
moments/mean/reduction_indices?
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/mean?
moments/StopGradientStopGradientmoments/mean:output:0*
T0**
_output_shapes
:2
moments/StopGradient?
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
moments/SquaredDifference?
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2$
"moments/variance/reduction_indices?
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/variance?
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze?
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze_1?
AssignMovingAvg/decayConst**
_class 
loc:@AssignMovingAvg/1222403*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg/decay?
AssignMovingAvg/ReadVariableOpReadVariableOpassignmovingavg_1222403*
_output_shapes
:*
dtype02 
AssignMovingAvg/ReadVariableOp?
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0**
_class 
loc:@AssignMovingAvg/1222403*
_output_shapes
:2
AssignMovingAvg/sub?
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0**
_class 
loc:@AssignMovingAvg/1222403*
_output_shapes
:2
AssignMovingAvg/mul?
#AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1222403AssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp**
_class 
loc:@AssignMovingAvg/1222403*
_output_shapes
 *
dtype02%
#AssignMovingAvg/AssignSubVariableOp?
AssignMovingAvg_1/decayConst*,
_class"
 loc:@AssignMovingAvg_1/1222409*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg_1/decay?
 AssignMovingAvg_1/ReadVariableOpReadVariableOpassignmovingavg_1_1222409*
_output_shapes
:*
dtype02"
 AssignMovingAvg_1/ReadVariableOp?
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1222409*
_output_shapes
:2
AssignMovingAvg_1/sub?
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1222409*
_output_shapes
:2
AssignMovingAvg_1/mul?
%AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1_1222409AssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*,
_class"
 loc:@AssignMovingAvg_1/1222409*
_output_shapes
 *
dtype02'
%AssignMovingAvg_1/AssignSubVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/mul_1{
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp?
batchnorm/subSub batchnorm/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0$^AssignMovingAvg/AssignSubVariableOp&^AssignMovingAvg_1/AssignSubVariableOp*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????::::2J
#AssignMovingAvg/AssignSubVariableOp#AssignMovingAvg/AssignSubVariableOp2N
%AssignMovingAvg_1/AssignSubVariableOp%AssignMovingAvg_1/AssignSubVariableOp:v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
?
T
(__inference_lambda_layer_call_fn_1226389
inputs_0
inputs_1
identity?
PartitionedCallPartitionedCallinputs_0inputs_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *L
fGRE
C__inference_lambda_layer_call_and_return_conditional_losses_12244842
PartitionedCalll
IdentityIdentityPartitionedCall:output:0*
T0*'
_output_shapes
:?????????
2

Identity"
identityIdentity:output:0*9
_input_shapes(
&:?????????
:?????????
:Q M
'
_output_shapes
:?????????

"
_user_specified_name
inputs/0:QM
'
_output_shapes
:?????????

"
_user_specified_name
inputs/1
?
?
.__inference_functional_3_layer_call_fn_1225055
input_1
input_2
input_3
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16

unknown_17

unknown_18

unknown_19

unknown_20

unknown_21

unknown_22

unknown_23

unknown_24

unknown_25

unknown_26

unknown_27

unknown_28

unknown_29

unknown_30
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinput_1input_2input_3unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19
unknown_20
unknown_21
unknown_22
unknown_23
unknown_24
unknown_25
unknown_26
unknown_27
unknown_28
unknown_29
unknown_30*.
Tin'
%2#*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????*B
_read_only_resource_inputs$
" 	
 !"*0
config_proto 

CPU

GPU2*0J 8? *R
fMRK
I__inference_functional_3_layer_call_and_return_conditional_losses_12249882
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????:??????????:?????????
::::::::::::::::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:] Y
4
_output_shapes"
 :??????????
!
_user_specified_name	input_1:]Y
4
_output_shapes"
 :??????????
!
_user_specified_name	input_2:PL
'
_output_shapes
:?????????

!
_user_specified_name	input_3
?	
?
E__inference_conv3d_1_layer_call_and_return_conditional_losses_1227104

inputs"
conv3d_readvariableop_resource#
biasadd_readvariableop_resource
identity??
Conv3D/ReadVariableOpReadVariableOpconv3d_readvariableop_resource**
_output_shapes
:*
dtype02
Conv3D/ReadVariableOp?
Conv3DConv3DinputsConv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
Conv3D?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddConv3D:output:0BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2	
BiasAdda
EluEluBiasAdd:output:0*
T0*3
_output_shapes!
:?????????2
Eluq
IdentityIdentityElu:activations:0*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*:
_input_shapes)
':?????????:::[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?+
?
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_1222708

inputs
assignmovingavg_1222683
assignmovingavg_1_1222689)
%batchnorm_mul_readvariableop_resource%
!batchnorm_readvariableop_resource
identity??#AssignMovingAvg/AssignSubVariableOp?%AssignMovingAvg_1/AssignSubVariableOp?
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2 
moments/mean/reduction_indices?
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/mean?
moments/StopGradientStopGradientmoments/mean:output:0*
T0**
_output_shapes
:2
moments/StopGradient?
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
moments/SquaredDifference?
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2$
"moments/variance/reduction_indices?
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/variance?
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze?
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze_1?
AssignMovingAvg/decayConst**
_class 
loc:@AssignMovingAvg/1222683*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg/decay?
AssignMovingAvg/ReadVariableOpReadVariableOpassignmovingavg_1222683*
_output_shapes
:*
dtype02 
AssignMovingAvg/ReadVariableOp?
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0**
_class 
loc:@AssignMovingAvg/1222683*
_output_shapes
:2
AssignMovingAvg/sub?
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0**
_class 
loc:@AssignMovingAvg/1222683*
_output_shapes
:2
AssignMovingAvg/mul?
#AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1222683AssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp**
_class 
loc:@AssignMovingAvg/1222683*
_output_shapes
 *
dtype02%
#AssignMovingAvg/AssignSubVariableOp?
AssignMovingAvg_1/decayConst*,
_class"
 loc:@AssignMovingAvg_1/1222689*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg_1/decay?
 AssignMovingAvg_1/ReadVariableOpReadVariableOpassignmovingavg_1_1222689*
_output_shapes
:*
dtype02"
 AssignMovingAvg_1/ReadVariableOp?
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1222689*
_output_shapes
:2
AssignMovingAvg_1/sub?
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1222689*
_output_shapes
:2
AssignMovingAvg_1/mul?
%AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1_1222689AssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*,
_class"
 loc:@AssignMovingAvg_1/1222689*
_output_shapes
 *
dtype02'
%AssignMovingAvg_1/AssignSubVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/mul_1{
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp?
batchnorm/subSub batchnorm/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0$^AssignMovingAvg/AssignSubVariableOp&^AssignMovingAvg_1/AssignSubVariableOp*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????::::2J
#AssignMovingAvg/AssignSubVariableOp#AssignMovingAvg/AssignSubVariableOp2N
%AssignMovingAvg_1/AssignSubVariableOp%AssignMovingAvg_1/AssignSubVariableOp:v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
?+
?
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_1222568

inputs
assignmovingavg_1222543
assignmovingavg_1_1222549)
%batchnorm_mul_readvariableop_resource%
!batchnorm_readvariableop_resource
identity??#AssignMovingAvg/AssignSubVariableOp?%AssignMovingAvg_1/AssignSubVariableOp?
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2 
moments/mean/reduction_indices?
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/mean?
moments/StopGradientStopGradientmoments/mean:output:0*
T0**
_output_shapes
:2
moments/StopGradient?
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
moments/SquaredDifference?
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2$
"moments/variance/reduction_indices?
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/variance?
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze?
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze_1?
AssignMovingAvg/decayConst**
_class 
loc:@AssignMovingAvg/1222543*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg/decay?
AssignMovingAvg/ReadVariableOpReadVariableOpassignmovingavg_1222543*
_output_shapes
:*
dtype02 
AssignMovingAvg/ReadVariableOp?
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0**
_class 
loc:@AssignMovingAvg/1222543*
_output_shapes
:2
AssignMovingAvg/sub?
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0**
_class 
loc:@AssignMovingAvg/1222543*
_output_shapes
:2
AssignMovingAvg/mul?
#AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1222543AssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp**
_class 
loc:@AssignMovingAvg/1222543*
_output_shapes
 *
dtype02%
#AssignMovingAvg/AssignSubVariableOp?
AssignMovingAvg_1/decayConst*,
_class"
 loc:@AssignMovingAvg_1/1222549*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg_1/decay?
 AssignMovingAvg_1/ReadVariableOpReadVariableOpassignmovingavg_1_1222549*
_output_shapes
:*
dtype02"
 AssignMovingAvg_1/ReadVariableOp?
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1222549*
_output_shapes
:2
AssignMovingAvg_1/sub?
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1222549*
_output_shapes
:2
AssignMovingAvg_1/mul?
%AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1_1222549AssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*,
_class"
 loc:@AssignMovingAvg_1/1222549*
_output_shapes
 *
dtype02'
%AssignMovingAvg_1/AssignSubVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/mul_1{
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp?
batchnorm/subSub batchnorm/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0$^AssignMovingAvg/AssignSubVariableOp&^AssignMovingAvg_1/AssignSubVariableOp*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????::::2J
#AssignMovingAvg/AssignSubVariableOp#AssignMovingAvg/AssignSubVariableOp2N
%AssignMovingAvg_1/AssignSubVariableOp%AssignMovingAvg_1/AssignSubVariableOp:v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
?
?
B__inference_dense_layer_call_and_return_conditional_losses_1223819

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity??
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	?d*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????d2
MatMul?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:d*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????d2	
BiasAddU
EluEluBiasAdd:output:0*
T0*'
_output_shapes
:?????????d2
Elu?
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	?d*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOp?
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?d2!
dense/kernel/Regularizer/Square?
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const?
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/Sum?
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2 
dense/kernel/Regularizer/mul/x?
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/mule
IdentityIdentityElu:activations:0*
T0*'
_output_shapes
:?????????d2

Identity"
identityIdentity:output:0*/
_input_shapes
:??????????:::P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs
?
?
7__inference_batch_normalization_2_layer_call_fn_1227379

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *N
_output_shapes<
::8????????????????????????????????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_12226012
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????::::22
StatefulPartitionedCallStatefulPartitionedCall:v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
?
o
C__inference_lambda_layer_call_and_return_conditional_losses_1226376
inputs_0
inputs_1
identityW
subSubinputs_0inputs_1*
T0*'
_output_shapes
:?????????
2
subL
AbsAbssub:z:0*
T0*'
_output_shapes
:?????????
2
Abs[
IdentityIdentityAbs:y:0*
T0*'
_output_shapes
:?????????
2

Identity"
identityIdentity:output:0*9
_input_shapes(
&:?????????
:?????????
:Q M
'
_output_shapes
:?????????

"
_user_specified_name
inputs/0:QM
'
_output_shapes
:?????????

"
_user_specified_name
inputs/1
?
?
P__inference_batch_normalization_layer_call_and_return_conditional_losses_1222849

inputs%
!batchnorm_readvariableop_resource)
%batchnorm_mul_readvariableop_resource'
#batchnorm_readvariableop_1_resource'
#batchnorm_readvariableop_2_resource
identity??
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/mul_1?
batchnorm/ReadVariableOp_1ReadVariableOp#batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_1?
batchnorm/mul_2Mul"batchnorm/ReadVariableOp_1:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOp_2ReadVariableOp#batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_2?
batchnorm/subSub"batchnorm/ReadVariableOp_2:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/add_1s
IdentityIdentitybatchnorm/add_1:z:0*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????:::::[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?
j
N__inference_average_pooling3d_layer_call_and_return_conditional_losses_1222758

inputs
identity?
	AvgPool3D	AvgPool3Dinputs*
T0*W
_output_shapesE
C:A?????????????????????????????????????????????*
ksize	
*
paddingVALID*
strides	
2
	AvgPool3D?
IdentityIdentityAvgPool3D:output:0*
T0*W
_output_shapesE
C:A?????????????????????????????????????????????2

Identity"
identityIdentity:output:0*V
_input_shapesE
C:A?????????????????????????????????????????????: {
W
_output_shapesE
C:A?????????????????????????????????????????????
 
_user_specified_nameinputs
?
?
D__inference_dense_1_layer_call_and_return_conditional_losses_1226879

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity??
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:d
*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????
2
MatMul?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:
*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????
2	
BiasAddU
EluEluBiasAdd:output:0*
T0*'
_output_shapes
:?????????
2
Elu?
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:d
*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp?
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d
2#
!dense_1/kernel/Regularizer/Square?
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const?
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/Sum?
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2"
 dense_1/kernel/Regularizer/mul/x?
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/mule
IdentityIdentityElu:activations:0*
T0*'
_output_shapes
:?????????
2

Identity"
identityIdentity:output:0*.
_input_shapes
:?????????d:::O K
'
_output_shapes
:?????????d
 
_user_specified_nameinputs
?	
?
E__inference_conv3d_2_layer_call_and_return_conditional_losses_1223014

inputs"
conv3d_readvariableop_resource#
biasadd_readvariableop_resource
identity??
Conv3D/ReadVariableOpReadVariableOpconv3d_readvariableop_resource**
_output_shapes
:*
dtype02
Conv3D/ReadVariableOp?
Conv3DConv3DinputsConv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
Conv3D?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddConv3D:output:0BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2	
BiasAdda
EluEluBiasAdd:output:0*
T0*3
_output_shapes!
:?????????2
Eluq
IdentityIdentityElu:activations:0*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*:
_input_shapes)
':?????????:::[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_1227619

inputs%
!batchnorm_readvariableop_resource)
%batchnorm_mul_readvariableop_resource'
#batchnorm_readvariableop_1_resource'
#batchnorm_readvariableop_2_resource
identity??
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/mul_1?
batchnorm/ReadVariableOp_1ReadVariableOp#batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_1?
batchnorm/mul_2Mul"batchnorm/ReadVariableOp_1:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOp_2ReadVariableOp#batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_2?
batchnorm/subSub"batchnorm/ReadVariableOp_2:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/add_1s
IdentityIdentitybatchnorm/add_1:z:0*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????:::::[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?

m
__inference_loss_fn_0_1226899;
7dense_kernel_regularizer_square_readvariableop_resource
identity??
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOp7dense_kernel_regularizer_square_readvariableop_resource*
_output_shapes
:	?d*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOp?
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?d2!
dense/kernel/Regularizer/Square?
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const?
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/Sum?
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2 
dense/kernel/Regularizer/mul/x?
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/mulc
IdentityIdentity dense/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: 2

Identity"
identityIdentity:output:0*
_input_shapes
:
?+
?
P__inference_batch_normalization_layer_call_and_return_conditional_losses_1226965

inputs
assignmovingavg_1226940
assignmovingavg_1_1226946)
%batchnorm_mul_readvariableop_resource%
!batchnorm_readvariableop_resource
identity??#AssignMovingAvg/AssignSubVariableOp?%AssignMovingAvg_1/AssignSubVariableOp?
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2 
moments/mean/reduction_indices?
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/mean?
moments/StopGradientStopGradientmoments/mean:output:0*
T0**
_output_shapes
:2
moments/StopGradient?
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
moments/SquaredDifference?
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2$
"moments/variance/reduction_indices?
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/variance?
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze?
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze_1?
AssignMovingAvg/decayConst**
_class 
loc:@AssignMovingAvg/1226940*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg/decay?
AssignMovingAvg/ReadVariableOpReadVariableOpassignmovingavg_1226940*
_output_shapes
:*
dtype02 
AssignMovingAvg/ReadVariableOp?
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0**
_class 
loc:@AssignMovingAvg/1226940*
_output_shapes
:2
AssignMovingAvg/sub?
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0**
_class 
loc:@AssignMovingAvg/1226940*
_output_shapes
:2
AssignMovingAvg/mul?
#AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1226940AssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp**
_class 
loc:@AssignMovingAvg/1226940*
_output_shapes
 *
dtype02%
#AssignMovingAvg/AssignSubVariableOp?
AssignMovingAvg_1/decayConst*,
_class"
 loc:@AssignMovingAvg_1/1226946*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg_1/decay?
 AssignMovingAvg_1/ReadVariableOpReadVariableOpassignmovingavg_1_1226946*
_output_shapes
:*
dtype02"
 AssignMovingAvg_1/ReadVariableOp?
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1226946*
_output_shapes
:2
AssignMovingAvg_1/sub?
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1226946*
_output_shapes
:2
AssignMovingAvg_1/mul?
%AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1_1226946AssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*,
_class"
 loc:@AssignMovingAvg_1/1226946*
_output_shapes
 *
dtype02'
%AssignMovingAvg_1/AssignSubVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/mul_1{
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp?
batchnorm/subSub batchnorm/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0$^AssignMovingAvg/AssignSubVariableOp&^AssignMovingAvg_1/AssignSubVariableOp*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????::::2J
#AssignMovingAvg/AssignSubVariableOp#AssignMovingAvg/AssignSubVariableOp2N
%AssignMovingAvg_1/AssignSubVariableOp%AssignMovingAvg_1/AssignSubVariableOp:v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
??
?&
 __inference__traced_save_1227989
file_prefix-
)savev2_dense_2_kernel_read_readvariableop+
'savev2_dense_2_bias_read_readvariableop(
$savev2_adam_iter_read_readvariableop	*
&savev2_adam_beta_1_read_readvariableop*
&savev2_adam_beta_2_read_readvariableop)
%savev2_adam_decay_read_readvariableop1
-savev2_adam_learning_rate_read_readvariableop,
(savev2_conv3d_kernel_read_readvariableop*
&savev2_conv3d_bias_read_readvariableop8
4savev2_batch_normalization_gamma_read_readvariableop7
3savev2_batch_normalization_beta_read_readvariableop.
*savev2_conv3d_1_kernel_read_readvariableop,
(savev2_conv3d_1_bias_read_readvariableop:
6savev2_batch_normalization_1_gamma_read_readvariableop9
5savev2_batch_normalization_1_beta_read_readvariableop.
*savev2_conv3d_2_kernel_read_readvariableop,
(savev2_conv3d_2_bias_read_readvariableop:
6savev2_batch_normalization_2_gamma_read_readvariableop9
5savev2_batch_normalization_2_beta_read_readvariableop.
*savev2_conv3d_3_kernel_read_readvariableop,
(savev2_conv3d_3_bias_read_readvariableop:
6savev2_batch_normalization_3_gamma_read_readvariableop9
5savev2_batch_normalization_3_beta_read_readvariableop,
(savev2_layer1_kernel_read_readvariableop*
&savev2_layer1_bias_read_readvariableop+
'savev2_dense_kernel_read_readvariableop)
%savev2_dense_bias_read_readvariableop-
)savev2_dense_1_kernel_read_readvariableop+
'savev2_dense_1_bias_read_readvariableop>
:savev2_batch_normalization_moving_mean_read_readvariableopB
>savev2_batch_normalization_moving_variance_read_readvariableop@
<savev2_batch_normalization_1_moving_mean_read_readvariableopD
@savev2_batch_normalization_1_moving_variance_read_readvariableop@
<savev2_batch_normalization_2_moving_mean_read_readvariableopD
@savev2_batch_normalization_2_moving_variance_read_readvariableop@
<savev2_batch_normalization_3_moving_mean_read_readvariableopD
@savev2_batch_normalization_3_moving_variance_read_readvariableop$
 savev2_total_read_readvariableop$
 savev2_count_read_readvariableop4
0savev2_adam_dense_2_kernel_m_read_readvariableop2
.savev2_adam_dense_2_bias_m_read_readvariableop3
/savev2_adam_conv3d_kernel_m_read_readvariableop1
-savev2_adam_conv3d_bias_m_read_readvariableop?
;savev2_adam_batch_normalization_gamma_m_read_readvariableop>
:savev2_adam_batch_normalization_beta_m_read_readvariableop5
1savev2_adam_conv3d_1_kernel_m_read_readvariableop3
/savev2_adam_conv3d_1_bias_m_read_readvariableopA
=savev2_adam_batch_normalization_1_gamma_m_read_readvariableop@
<savev2_adam_batch_normalization_1_beta_m_read_readvariableop5
1savev2_adam_conv3d_2_kernel_m_read_readvariableop3
/savev2_adam_conv3d_2_bias_m_read_readvariableopA
=savev2_adam_batch_normalization_2_gamma_m_read_readvariableop@
<savev2_adam_batch_normalization_2_beta_m_read_readvariableop5
1savev2_adam_conv3d_3_kernel_m_read_readvariableop3
/savev2_adam_conv3d_3_bias_m_read_readvariableopA
=savev2_adam_batch_normalization_3_gamma_m_read_readvariableop@
<savev2_adam_batch_normalization_3_beta_m_read_readvariableop3
/savev2_adam_layer1_kernel_m_read_readvariableop1
-savev2_adam_layer1_bias_m_read_readvariableop2
.savev2_adam_dense_kernel_m_read_readvariableop0
,savev2_adam_dense_bias_m_read_readvariableop4
0savev2_adam_dense_1_kernel_m_read_readvariableop2
.savev2_adam_dense_1_bias_m_read_readvariableop4
0savev2_adam_dense_2_kernel_v_read_readvariableop2
.savev2_adam_dense_2_bias_v_read_readvariableop3
/savev2_adam_conv3d_kernel_v_read_readvariableop1
-savev2_adam_conv3d_bias_v_read_readvariableop?
;savev2_adam_batch_normalization_gamma_v_read_readvariableop>
:savev2_adam_batch_normalization_beta_v_read_readvariableop5
1savev2_adam_conv3d_1_kernel_v_read_readvariableop3
/savev2_adam_conv3d_1_bias_v_read_readvariableopA
=savev2_adam_batch_normalization_1_gamma_v_read_readvariableop@
<savev2_adam_batch_normalization_1_beta_v_read_readvariableop5
1savev2_adam_conv3d_2_kernel_v_read_readvariableop3
/savev2_adam_conv3d_2_bias_v_read_readvariableopA
=savev2_adam_batch_normalization_2_gamma_v_read_readvariableop@
<savev2_adam_batch_normalization_2_beta_v_read_readvariableop5
1savev2_adam_conv3d_3_kernel_v_read_readvariableop3
/savev2_adam_conv3d_3_bias_v_read_readvariableopA
=savev2_adam_batch_normalization_3_gamma_v_read_readvariableop@
<savev2_adam_batch_normalization_3_beta_v_read_readvariableop3
/savev2_adam_layer1_kernel_v_read_readvariableop1
-savev2_adam_layer1_bias_v_read_readvariableop2
.savev2_adam_dense_kernel_v_read_readvariableop0
,savev2_adam_dense_bias_v_read_readvariableop4
0savev2_adam_dense_1_kernel_v_read_readvariableop2
.savev2_adam_dense_1_bias_v_read_readvariableop
savev2_const

identity_1??MergeV2Checkpoints?
StaticRegexFullMatchStaticRegexFullMatchfile_prefix"/device:CPU:**
_output_shapes
: *
pattern
^s3://.*2
StaticRegexFullMatchc
ConstConst"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B.part2
Const?
Const_1Const"/device:CPU:**
_output_shapes
: *
dtype0*<
value3B1 B+_temp_bd49cba5093c41b484e85447e6dfc7da/part2	
Const_1?
SelectSelectStaticRegexFullMatch:output:0Const:output:0Const_1:output:0"/device:CPU:**
T0*
_output_shapes
: 2
Selectt

StringJoin
StringJoinfile_prefixSelect:output:0"/device:CPU:**
N*
_output_shapes
: 2

StringJoinZ

num_shardsConst*
_output_shapes
: *
dtype0*
value	B :2

num_shards
ShardedFilename/shardConst"/device:CPU:0*
_output_shapes
: *
dtype0*
value	B : 2
ShardedFilename/shard?
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: 2
ShardedFilename?-
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:X*
dtype0*?,
value?,B?,XB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/0/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/1/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/2/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/3/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/4/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/5/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/6/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/7/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/8/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/9/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/10/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/11/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/12/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/13/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/14/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/15/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/16/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/17/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/18/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/19/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/20/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/21/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB'variables/10/.ATTRIBUTES/VARIABLE_VALUEB'variables/11/.ATTRIBUTES/VARIABLE_VALUEB'variables/16/.ATTRIBUTES/VARIABLE_VALUEB'variables/17/.ATTRIBUTES/VARIABLE_VALUEB'variables/22/.ATTRIBUTES/VARIABLE_VALUEB'variables/23/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/0/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/1/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/4/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/5/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/7/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/8/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/9/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/10/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/11/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/12/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/13/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/14/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/15/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/16/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/17/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/18/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/19/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/20/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/21/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/0/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/1/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/4/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/5/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/7/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/8/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/9/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/10/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/11/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/12/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/13/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/14/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/15/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/16/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/17/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/18/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/19/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/20/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/21/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH2
SaveV2/tensor_names?
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:X*
dtype0*?
value?B?XB B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 2
SaveV2/shape_and_slices?$
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0)savev2_dense_2_kernel_read_readvariableop'savev2_dense_2_bias_read_readvariableop$savev2_adam_iter_read_readvariableop&savev2_adam_beta_1_read_readvariableop&savev2_adam_beta_2_read_readvariableop%savev2_adam_decay_read_readvariableop-savev2_adam_learning_rate_read_readvariableop(savev2_conv3d_kernel_read_readvariableop&savev2_conv3d_bias_read_readvariableop4savev2_batch_normalization_gamma_read_readvariableop3savev2_batch_normalization_beta_read_readvariableop*savev2_conv3d_1_kernel_read_readvariableop(savev2_conv3d_1_bias_read_readvariableop6savev2_batch_normalization_1_gamma_read_readvariableop5savev2_batch_normalization_1_beta_read_readvariableop*savev2_conv3d_2_kernel_read_readvariableop(savev2_conv3d_2_bias_read_readvariableop6savev2_batch_normalization_2_gamma_read_readvariableop5savev2_batch_normalization_2_beta_read_readvariableop*savev2_conv3d_3_kernel_read_readvariableop(savev2_conv3d_3_bias_read_readvariableop6savev2_batch_normalization_3_gamma_read_readvariableop5savev2_batch_normalization_3_beta_read_readvariableop(savev2_layer1_kernel_read_readvariableop&savev2_layer1_bias_read_readvariableop'savev2_dense_kernel_read_readvariableop%savev2_dense_bias_read_readvariableop)savev2_dense_1_kernel_read_readvariableop'savev2_dense_1_bias_read_readvariableop:savev2_batch_normalization_moving_mean_read_readvariableop>savev2_batch_normalization_moving_variance_read_readvariableop<savev2_batch_normalization_1_moving_mean_read_readvariableop@savev2_batch_normalization_1_moving_variance_read_readvariableop<savev2_batch_normalization_2_moving_mean_read_readvariableop@savev2_batch_normalization_2_moving_variance_read_readvariableop<savev2_batch_normalization_3_moving_mean_read_readvariableop@savev2_batch_normalization_3_moving_variance_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableop0savev2_adam_dense_2_kernel_m_read_readvariableop.savev2_adam_dense_2_bias_m_read_readvariableop/savev2_adam_conv3d_kernel_m_read_readvariableop-savev2_adam_conv3d_bias_m_read_readvariableop;savev2_adam_batch_normalization_gamma_m_read_readvariableop:savev2_adam_batch_normalization_beta_m_read_readvariableop1savev2_adam_conv3d_1_kernel_m_read_readvariableop/savev2_adam_conv3d_1_bias_m_read_readvariableop=savev2_adam_batch_normalization_1_gamma_m_read_readvariableop<savev2_adam_batch_normalization_1_beta_m_read_readvariableop1savev2_adam_conv3d_2_kernel_m_read_readvariableop/savev2_adam_conv3d_2_bias_m_read_readvariableop=savev2_adam_batch_normalization_2_gamma_m_read_readvariableop<savev2_adam_batch_normalization_2_beta_m_read_readvariableop1savev2_adam_conv3d_3_kernel_m_read_readvariableop/savev2_adam_conv3d_3_bias_m_read_readvariableop=savev2_adam_batch_normalization_3_gamma_m_read_readvariableop<savev2_adam_batch_normalization_3_beta_m_read_readvariableop/savev2_adam_layer1_kernel_m_read_readvariableop-savev2_adam_layer1_bias_m_read_readvariableop.savev2_adam_dense_kernel_m_read_readvariableop,savev2_adam_dense_bias_m_read_readvariableop0savev2_adam_dense_1_kernel_m_read_readvariableop.savev2_adam_dense_1_bias_m_read_readvariableop0savev2_adam_dense_2_kernel_v_read_readvariableop.savev2_adam_dense_2_bias_v_read_readvariableop/savev2_adam_conv3d_kernel_v_read_readvariableop-savev2_adam_conv3d_bias_v_read_readvariableop;savev2_adam_batch_normalization_gamma_v_read_readvariableop:savev2_adam_batch_normalization_beta_v_read_readvariableop1savev2_adam_conv3d_1_kernel_v_read_readvariableop/savev2_adam_conv3d_1_bias_v_read_readvariableop=savev2_adam_batch_normalization_1_gamma_v_read_readvariableop<savev2_adam_batch_normalization_1_beta_v_read_readvariableop1savev2_adam_conv3d_2_kernel_v_read_readvariableop/savev2_adam_conv3d_2_bias_v_read_readvariableop=savev2_adam_batch_normalization_2_gamma_v_read_readvariableop<savev2_adam_batch_normalization_2_beta_v_read_readvariableop1savev2_adam_conv3d_3_kernel_v_read_readvariableop/savev2_adam_conv3d_3_bias_v_read_readvariableop=savev2_adam_batch_normalization_3_gamma_v_read_readvariableop<savev2_adam_batch_normalization_3_beta_v_read_readvariableop/savev2_adam_layer1_kernel_v_read_readvariableop-savev2_adam_layer1_bias_v_read_readvariableop.savev2_adam_dense_kernel_v_read_readvariableop,savev2_adam_dense_bias_v_read_readvariableop0savev2_adam_dense_1_kernel_v_read_readvariableop.savev2_adam_dense_1_bias_v_read_readvariableopsavev2_const"/device:CPU:0*
_output_shapes
 *f
dtypes\
Z2X	2
SaveV2?
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0^SaveV2"/device:CPU:0*
N*
T0*
_output_shapes
:2(
&MergeV2Checkpoints/checkpoint_prefixes?
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix"/device:CPU:0*
_output_shapes
 2
MergeV2Checkpointsr
IdentityIdentityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: 2

Identitym

Identity_1IdentityIdentity:output:0^MergeV2Checkpoints*
T0*
_output_shapes
: 2

Identity_1"!

identity_1Identity_1:output:0*?
_input_shapes?
?: ::: : : : : :?::::::::::::::::
?
?:?:	?d:d:d
:
::::::::: : :::?::::::::::::::::
?
?:?:	?d:d:d
:
:::?::::::::::::::::
?
?:?:	?d:d:d
:
: 2(
MergeV2CheckpointsMergeV2Checkpoints:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix:$ 

_output_shapes

:: 

_output_shapes
::

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :1-
+
_output_shapes
:?: 	

_output_shapes
:: 


_output_shapes
:: 

_output_shapes
::0,
*
_output_shapes
:: 

_output_shapes
:: 

_output_shapes
:: 

_output_shapes
::0,
*
_output_shapes
:: 

_output_shapes
:: 

_output_shapes
:: 

_output_shapes
::0,
*
_output_shapes
:: 

_output_shapes
:: 

_output_shapes
:: 

_output_shapes
::&"
 
_output_shapes
:
?
?:!

_output_shapes	
:?:%!

_output_shapes
:	?d: 

_output_shapes
:d:$ 

_output_shapes

:d
: 

_output_shapes
:
: 

_output_shapes
:: 

_output_shapes
::  

_output_shapes
:: !

_output_shapes
:: "

_output_shapes
:: #

_output_shapes
:: $

_output_shapes
:: %

_output_shapes
::&

_output_shapes
: :'

_output_shapes
: :$( 

_output_shapes

:: )

_output_shapes
::1*-
+
_output_shapes
:?: +

_output_shapes
:: ,

_output_shapes
:: -

_output_shapes
::0.,
*
_output_shapes
:: /

_output_shapes
:: 0

_output_shapes
:: 1

_output_shapes
::02,
*
_output_shapes
:: 3

_output_shapes
:: 4

_output_shapes
:: 5

_output_shapes
::06,
*
_output_shapes
:: 7

_output_shapes
:: 8

_output_shapes
:: 9

_output_shapes
::&:"
 
_output_shapes
:
?
?:!;

_output_shapes	
:?:%<!

_output_shapes
:	?d: =

_output_shapes
:d:$> 

_output_shapes

:d
: ?

_output_shapes
:
:$@ 

_output_shapes

:: A

_output_shapes
::1B-
+
_output_shapes
:?: C

_output_shapes
:: D

_output_shapes
:: E

_output_shapes
::0F,
*
_output_shapes
:: G

_output_shapes
:: H

_output_shapes
:: I

_output_shapes
::0J,
*
_output_shapes
:: K

_output_shapes
:: L

_output_shapes
:: M

_output_shapes
::0N,
*
_output_shapes
:: O

_output_shapes
:: P

_output_shapes
:: Q

_output_shapes
::&R"
 
_output_shapes
:
?
?:!S

_output_shapes	
:?:%T!

_output_shapes
:	?d: U

_output_shapes
:d:$V 

_output_shapes

:d
: W

_output_shapes
:
:X

_output_shapes
: 
?
E
)__inference_dropout_layer_call_fn_1227683

inputs
identity?
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????
* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dropout_layer_call_and_return_conditional_losses_12232712
PartitionedCallm
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:??????????
2

Identity"
identityIdentity:output:0*'
_input_shapes
:??????????
:P L
(
_output_shapes
:??????????

 
_user_specified_nameinputs
?

*__inference_conv3d_1_layer_call_fn_1227113

inputs
unknown
	unknown_0
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv3d_1_layer_call_and_return_conditional_losses_12228962
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*:
_input_shapes)
':?????????::22
StatefulPartitionedCallStatefulPartitionedCall:[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?

*__inference_conv3d_2_layer_call_fn_1227297

inputs
unknown
	unknown_0
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv3d_2_layer_call_and_return_conditional_losses_12230142
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*:
_input_shapes)
':?????????::22
StatefulPartitionedCallStatefulPartitionedCall:[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?
?
5__inference_batch_normalization_layer_call_fn_1227011

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *N
_output_shapes<
::8????????????????????????????????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *Y
fTRR
P__inference_batch_normalization_layer_call_and_return_conditional_losses_12223212
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????::::22
StatefulPartitionedCallStatefulPartitionedCall:v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
?
T
(__inference_lambda_layer_call_fn_1226395
inputs_0
inputs_1
identity?
PartitionedCallPartitionedCallinputs_0inputs_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *L
fGRE
C__inference_lambda_layer_call_and_return_conditional_losses_12244912
PartitionedCalll
IdentityIdentityPartitionedCall:output:0*
T0*'
_output_shapes
:?????????
2

Identity"
identityIdentity:output:0*9
_input_shapes(
&:?????????
:?????????
:Q M
'
_output_shapes
:?????????

"
_user_specified_name
inputs/0:QM
'
_output_shapes
:?????????

"
_user_specified_name
inputs/1
?	
?
E__inference_conv3d_3_layer_call_and_return_conditional_losses_1227472

inputs"
conv3d_readvariableop_resource#
biasadd_readvariableop_resource
identity??
Conv3D/ReadVariableOpReadVariableOpconv3d_readvariableop_resource**
_output_shapes
:*
dtype02
Conv3D/ReadVariableOp?
Conv3DConv3DinputsConv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
Conv3D?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddConv3D:output:0BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2	
BiasAdda
EluEluBiasAdd:output:0*
T0*3
_output_shapes!
:?????????2
Eluq
IdentityIdentityElu:activations:0*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*:
_input_shapes)
':?????????:::[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
??
?
I__inference_functional_3_layer_call_and_return_conditional_losses_1225750
inputs_0
inputs_1
inputs_2A
=sequential_functional_1_conv3d_conv3d_readvariableop_resourceB
>sequential_functional_1_conv3d_biasadd_readvariableop_resourceQ
Msequential_functional_1_batch_normalization_batchnorm_readvariableop_resourceU
Qsequential_functional_1_batch_normalization_batchnorm_mul_readvariableop_resourceS
Osequential_functional_1_batch_normalization_batchnorm_readvariableop_1_resourceS
Osequential_functional_1_batch_normalization_batchnorm_readvariableop_2_resourceC
?sequential_functional_1_conv3d_1_conv3d_readvariableop_resourceD
@sequential_functional_1_conv3d_1_biasadd_readvariableop_resourceS
Osequential_functional_1_batch_normalization_1_batchnorm_readvariableop_resourceW
Ssequential_functional_1_batch_normalization_1_batchnorm_mul_readvariableop_resourceU
Qsequential_functional_1_batch_normalization_1_batchnorm_readvariableop_1_resourceU
Qsequential_functional_1_batch_normalization_1_batchnorm_readvariableop_2_resourceC
?sequential_functional_1_conv3d_2_conv3d_readvariableop_resourceD
@sequential_functional_1_conv3d_2_biasadd_readvariableop_resourceS
Osequential_functional_1_batch_normalization_2_batchnorm_readvariableop_resourceW
Ssequential_functional_1_batch_normalization_2_batchnorm_mul_readvariableop_resourceU
Qsequential_functional_1_batch_normalization_2_batchnorm_readvariableop_1_resourceU
Qsequential_functional_1_batch_normalization_2_batchnorm_readvariableop_2_resourceC
?sequential_functional_1_conv3d_3_conv3d_readvariableop_resourceD
@sequential_functional_1_conv3d_3_biasadd_readvariableop_resourceS
Osequential_functional_1_batch_normalization_3_batchnorm_readvariableop_resourceW
Ssequential_functional_1_batch_normalization_3_batchnorm_mul_readvariableop_resourceU
Qsequential_functional_1_batch_normalization_3_batchnorm_readvariableop_1_resourceU
Qsequential_functional_1_batch_normalization_3_batchnorm_readvariableop_2_resourceA
=sequential_functional_1_layer1_matmul_readvariableop_resourceB
>sequential_functional_1_layer1_biasadd_readvariableop_resource3
/sequential_dense_matmul_readvariableop_resource4
0sequential_dense_biasadd_readvariableop_resource5
1sequential_dense_1_matmul_readvariableop_resource6
2sequential_dense_1_biasadd_readvariableop_resource*
&dense_2_matmul_readvariableop_resource+
'dense_2_biasadd_readvariableop_resource
identity??
4sequential/functional_1/conv3d/Conv3D/ReadVariableOpReadVariableOp=sequential_functional_1_conv3d_conv3d_readvariableop_resource*+
_output_shapes
:?*
dtype026
4sequential/functional_1/conv3d/Conv3D/ReadVariableOp?
%sequential/functional_1/conv3d/Conv3DConv3Dinputs_0<sequential/functional_1/conv3d/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2'
%sequential/functional_1/conv3d/Conv3D?
5sequential/functional_1/conv3d/BiasAdd/ReadVariableOpReadVariableOp>sequential_functional_1_conv3d_biasadd_readvariableop_resource*
_output_shapes
:*
dtype027
5sequential/functional_1/conv3d/BiasAdd/ReadVariableOp?
&sequential/functional_1/conv3d/BiasAddBiasAdd.sequential/functional_1/conv3d/Conv3D:output:0=sequential/functional_1/conv3d/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2(
&sequential/functional_1/conv3d/BiasAdd?
Dsequential/functional_1/batch_normalization/batchnorm/ReadVariableOpReadVariableOpMsequential_functional_1_batch_normalization_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02F
Dsequential/functional_1/batch_normalization/batchnorm/ReadVariableOp?
;sequential/functional_1/batch_normalization/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2=
;sequential/functional_1/batch_normalization/batchnorm/add/y?
9sequential/functional_1/batch_normalization/batchnorm/addAddV2Lsequential/functional_1/batch_normalization/batchnorm/ReadVariableOp:value:0Dsequential/functional_1/batch_normalization/batchnorm/add/y:output:0*
T0*
_output_shapes
:2;
9sequential/functional_1/batch_normalization/batchnorm/add?
;sequential/functional_1/batch_normalization/batchnorm/RsqrtRsqrt=sequential/functional_1/batch_normalization/batchnorm/add:z:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization/batchnorm/Rsqrt?
Hsequential/functional_1/batch_normalization/batchnorm/mul/ReadVariableOpReadVariableOpQsequential_functional_1_batch_normalization_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02J
Hsequential/functional_1/batch_normalization/batchnorm/mul/ReadVariableOp?
9sequential/functional_1/batch_normalization/batchnorm/mulMul?sequential/functional_1/batch_normalization/batchnorm/Rsqrt:y:0Psequential/functional_1/batch_normalization/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2;
9sequential/functional_1/batch_normalization/batchnorm/mul?
;sequential/functional_1/batch_normalization/batchnorm/mul_1Mul/sequential/functional_1/conv3d/BiasAdd:output:0=sequential/functional_1/batch_normalization/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2=
;sequential/functional_1/batch_normalization/batchnorm/mul_1?
Fsequential/functional_1/batch_normalization/batchnorm/ReadVariableOp_1ReadVariableOpOsequential_functional_1_batch_normalization_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02H
Fsequential/functional_1/batch_normalization/batchnorm/ReadVariableOp_1?
;sequential/functional_1/batch_normalization/batchnorm/mul_2MulNsequential/functional_1/batch_normalization/batchnorm/ReadVariableOp_1:value:0=sequential/functional_1/batch_normalization/batchnorm/mul:z:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization/batchnorm/mul_2?
Fsequential/functional_1/batch_normalization/batchnorm/ReadVariableOp_2ReadVariableOpOsequential_functional_1_batch_normalization_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02H
Fsequential/functional_1/batch_normalization/batchnorm/ReadVariableOp_2?
9sequential/functional_1/batch_normalization/batchnorm/subSubNsequential/functional_1/batch_normalization/batchnorm/ReadVariableOp_2:value:0?sequential/functional_1/batch_normalization/batchnorm/mul_2:z:0*
T0*
_output_shapes
:2;
9sequential/functional_1/batch_normalization/batchnorm/sub?
;sequential/functional_1/batch_normalization/batchnorm/add_1AddV2?sequential/functional_1/batch_normalization/batchnorm/mul_1:z:0=sequential/functional_1/batch_normalization/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2=
;sequential/functional_1/batch_normalization/batchnorm/add_1?
6sequential/functional_1/conv3d_1/Conv3D/ReadVariableOpReadVariableOp?sequential_functional_1_conv3d_1_conv3d_readvariableop_resource**
_output_shapes
:*
dtype028
6sequential/functional_1/conv3d_1/Conv3D/ReadVariableOp?
'sequential/functional_1/conv3d_1/Conv3DConv3D?sequential/functional_1/batch_normalization/batchnorm/add_1:z:0>sequential/functional_1/conv3d_1/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2)
'sequential/functional_1/conv3d_1/Conv3D?
7sequential/functional_1/conv3d_1/BiasAdd/ReadVariableOpReadVariableOp@sequential_functional_1_conv3d_1_biasadd_readvariableop_resource*
_output_shapes
:*
dtype029
7sequential/functional_1/conv3d_1/BiasAdd/ReadVariableOp?
(sequential/functional_1/conv3d_1/BiasAddBiasAdd0sequential/functional_1/conv3d_1/Conv3D:output:0?sequential/functional_1/conv3d_1/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2*
(sequential/functional_1/conv3d_1/BiasAdd?
$sequential/functional_1/conv3d_1/EluElu1sequential/functional_1/conv3d_1/BiasAdd:output:0*
T0*3
_output_shapes!
:?????????2&
$sequential/functional_1/conv3d_1/Elu?
Fsequential/functional_1/batch_normalization_1/batchnorm/ReadVariableOpReadVariableOpOsequential_functional_1_batch_normalization_1_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02H
Fsequential/functional_1/batch_normalization_1/batchnorm/ReadVariableOp?
=sequential/functional_1/batch_normalization_1/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2?
=sequential/functional_1/batch_normalization_1/batchnorm/add/y?
;sequential/functional_1/batch_normalization_1/batchnorm/addAddV2Nsequential/functional_1/batch_normalization_1/batchnorm/ReadVariableOp:value:0Fsequential/functional_1/batch_normalization_1/batchnorm/add/y:output:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization_1/batchnorm/add?
=sequential/functional_1/batch_normalization_1/batchnorm/RsqrtRsqrt?sequential/functional_1/batch_normalization_1/batchnorm/add:z:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_1/batchnorm/Rsqrt?
Jsequential/functional_1/batch_normalization_1/batchnorm/mul/ReadVariableOpReadVariableOpSsequential_functional_1_batch_normalization_1_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02L
Jsequential/functional_1/batch_normalization_1/batchnorm/mul/ReadVariableOp?
;sequential/functional_1/batch_normalization_1/batchnorm/mulMulAsequential/functional_1/batch_normalization_1/batchnorm/Rsqrt:y:0Rsequential/functional_1/batch_normalization_1/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization_1/batchnorm/mul?
=sequential/functional_1/batch_normalization_1/batchnorm/mul_1Mul2sequential/functional_1/conv3d_1/Elu:activations:0?sequential/functional_1/batch_normalization_1/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2?
=sequential/functional_1/batch_normalization_1/batchnorm/mul_1?
Hsequential/functional_1/batch_normalization_1/batchnorm/ReadVariableOp_1ReadVariableOpQsequential_functional_1_batch_normalization_1_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02J
Hsequential/functional_1/batch_normalization_1/batchnorm/ReadVariableOp_1?
=sequential/functional_1/batch_normalization_1/batchnorm/mul_2MulPsequential/functional_1/batch_normalization_1/batchnorm/ReadVariableOp_1:value:0?sequential/functional_1/batch_normalization_1/batchnorm/mul:z:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_1/batchnorm/mul_2?
Hsequential/functional_1/batch_normalization_1/batchnorm/ReadVariableOp_2ReadVariableOpQsequential_functional_1_batch_normalization_1_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02J
Hsequential/functional_1/batch_normalization_1/batchnorm/ReadVariableOp_2?
;sequential/functional_1/batch_normalization_1/batchnorm/subSubPsequential/functional_1/batch_normalization_1/batchnorm/ReadVariableOp_2:value:0Asequential/functional_1/batch_normalization_1/batchnorm/mul_2:z:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization_1/batchnorm/sub?
=sequential/functional_1/batch_normalization_1/batchnorm/add_1AddV2Asequential/functional_1/batch_normalization_1/batchnorm/mul_1:z:0?sequential/functional_1/batch_normalization_1/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2?
=sequential/functional_1/batch_normalization_1/batchnorm/add_1?
6sequential/functional_1/conv3d_2/Conv3D/ReadVariableOpReadVariableOp?sequential_functional_1_conv3d_2_conv3d_readvariableop_resource**
_output_shapes
:*
dtype028
6sequential/functional_1/conv3d_2/Conv3D/ReadVariableOp?
'sequential/functional_1/conv3d_2/Conv3DConv3DAsequential/functional_1/batch_normalization_1/batchnorm/add_1:z:0>sequential/functional_1/conv3d_2/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2)
'sequential/functional_1/conv3d_2/Conv3D?
7sequential/functional_1/conv3d_2/BiasAdd/ReadVariableOpReadVariableOp@sequential_functional_1_conv3d_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype029
7sequential/functional_1/conv3d_2/BiasAdd/ReadVariableOp?
(sequential/functional_1/conv3d_2/BiasAddBiasAdd0sequential/functional_1/conv3d_2/Conv3D:output:0?sequential/functional_1/conv3d_2/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2*
(sequential/functional_1/conv3d_2/BiasAdd?
$sequential/functional_1/conv3d_2/EluElu1sequential/functional_1/conv3d_2/BiasAdd:output:0*
T0*3
_output_shapes!
:?????????2&
$sequential/functional_1/conv3d_2/Elu?
Fsequential/functional_1/batch_normalization_2/batchnorm/ReadVariableOpReadVariableOpOsequential_functional_1_batch_normalization_2_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02H
Fsequential/functional_1/batch_normalization_2/batchnorm/ReadVariableOp?
=sequential/functional_1/batch_normalization_2/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2?
=sequential/functional_1/batch_normalization_2/batchnorm/add/y?
;sequential/functional_1/batch_normalization_2/batchnorm/addAddV2Nsequential/functional_1/batch_normalization_2/batchnorm/ReadVariableOp:value:0Fsequential/functional_1/batch_normalization_2/batchnorm/add/y:output:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization_2/batchnorm/add?
=sequential/functional_1/batch_normalization_2/batchnorm/RsqrtRsqrt?sequential/functional_1/batch_normalization_2/batchnorm/add:z:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_2/batchnorm/Rsqrt?
Jsequential/functional_1/batch_normalization_2/batchnorm/mul/ReadVariableOpReadVariableOpSsequential_functional_1_batch_normalization_2_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02L
Jsequential/functional_1/batch_normalization_2/batchnorm/mul/ReadVariableOp?
;sequential/functional_1/batch_normalization_2/batchnorm/mulMulAsequential/functional_1/batch_normalization_2/batchnorm/Rsqrt:y:0Rsequential/functional_1/batch_normalization_2/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization_2/batchnorm/mul?
=sequential/functional_1/batch_normalization_2/batchnorm/mul_1Mul2sequential/functional_1/conv3d_2/Elu:activations:0?sequential/functional_1/batch_normalization_2/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2?
=sequential/functional_1/batch_normalization_2/batchnorm/mul_1?
Hsequential/functional_1/batch_normalization_2/batchnorm/ReadVariableOp_1ReadVariableOpQsequential_functional_1_batch_normalization_2_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02J
Hsequential/functional_1/batch_normalization_2/batchnorm/ReadVariableOp_1?
=sequential/functional_1/batch_normalization_2/batchnorm/mul_2MulPsequential/functional_1/batch_normalization_2/batchnorm/ReadVariableOp_1:value:0?sequential/functional_1/batch_normalization_2/batchnorm/mul:z:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_2/batchnorm/mul_2?
Hsequential/functional_1/batch_normalization_2/batchnorm/ReadVariableOp_2ReadVariableOpQsequential_functional_1_batch_normalization_2_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02J
Hsequential/functional_1/batch_normalization_2/batchnorm/ReadVariableOp_2?
;sequential/functional_1/batch_normalization_2/batchnorm/subSubPsequential/functional_1/batch_normalization_2/batchnorm/ReadVariableOp_2:value:0Asequential/functional_1/batch_normalization_2/batchnorm/mul_2:z:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization_2/batchnorm/sub?
=sequential/functional_1/batch_normalization_2/batchnorm/add_1AddV2Asequential/functional_1/batch_normalization_2/batchnorm/mul_1:z:0?sequential/functional_1/batch_normalization_2/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2?
=sequential/functional_1/batch_normalization_2/batchnorm/add_1?
6sequential/functional_1/conv3d_3/Conv3D/ReadVariableOpReadVariableOp?sequential_functional_1_conv3d_3_conv3d_readvariableop_resource**
_output_shapes
:*
dtype028
6sequential/functional_1/conv3d_3/Conv3D/ReadVariableOp?
'sequential/functional_1/conv3d_3/Conv3DConv3DAsequential/functional_1/batch_normalization_2/batchnorm/add_1:z:0>sequential/functional_1/conv3d_3/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2)
'sequential/functional_1/conv3d_3/Conv3D?
7sequential/functional_1/conv3d_3/BiasAdd/ReadVariableOpReadVariableOp@sequential_functional_1_conv3d_3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype029
7sequential/functional_1/conv3d_3/BiasAdd/ReadVariableOp?
(sequential/functional_1/conv3d_3/BiasAddBiasAdd0sequential/functional_1/conv3d_3/Conv3D:output:0?sequential/functional_1/conv3d_3/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2*
(sequential/functional_1/conv3d_3/BiasAdd?
$sequential/functional_1/conv3d_3/EluElu1sequential/functional_1/conv3d_3/BiasAdd:output:0*
T0*3
_output_shapes!
:?????????2&
$sequential/functional_1/conv3d_3/Elu?
Fsequential/functional_1/batch_normalization_3/batchnorm/ReadVariableOpReadVariableOpOsequential_functional_1_batch_normalization_3_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02H
Fsequential/functional_1/batch_normalization_3/batchnorm/ReadVariableOp?
=sequential/functional_1/batch_normalization_3/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2?
=sequential/functional_1/batch_normalization_3/batchnorm/add/y?
;sequential/functional_1/batch_normalization_3/batchnorm/addAddV2Nsequential/functional_1/batch_normalization_3/batchnorm/ReadVariableOp:value:0Fsequential/functional_1/batch_normalization_3/batchnorm/add/y:output:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization_3/batchnorm/add?
=sequential/functional_1/batch_normalization_3/batchnorm/RsqrtRsqrt?sequential/functional_1/batch_normalization_3/batchnorm/add:z:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_3/batchnorm/Rsqrt?
Jsequential/functional_1/batch_normalization_3/batchnorm/mul/ReadVariableOpReadVariableOpSsequential_functional_1_batch_normalization_3_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02L
Jsequential/functional_1/batch_normalization_3/batchnorm/mul/ReadVariableOp?
;sequential/functional_1/batch_normalization_3/batchnorm/mulMulAsequential/functional_1/batch_normalization_3/batchnorm/Rsqrt:y:0Rsequential/functional_1/batch_normalization_3/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization_3/batchnorm/mul?
=sequential/functional_1/batch_normalization_3/batchnorm/mul_1Mul2sequential/functional_1/conv3d_3/Elu:activations:0?sequential/functional_1/batch_normalization_3/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2?
=sequential/functional_1/batch_normalization_3/batchnorm/mul_1?
Hsequential/functional_1/batch_normalization_3/batchnorm/ReadVariableOp_1ReadVariableOpQsequential_functional_1_batch_normalization_3_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02J
Hsequential/functional_1/batch_normalization_3/batchnorm/ReadVariableOp_1?
=sequential/functional_1/batch_normalization_3/batchnorm/mul_2MulPsequential/functional_1/batch_normalization_3/batchnorm/ReadVariableOp_1:value:0?sequential/functional_1/batch_normalization_3/batchnorm/mul:z:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_3/batchnorm/mul_2?
Hsequential/functional_1/batch_normalization_3/batchnorm/ReadVariableOp_2ReadVariableOpQsequential_functional_1_batch_normalization_3_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02J
Hsequential/functional_1/batch_normalization_3/batchnorm/ReadVariableOp_2?
;sequential/functional_1/batch_normalization_3/batchnorm/subSubPsequential/functional_1/batch_normalization_3/batchnorm/ReadVariableOp_2:value:0Asequential/functional_1/batch_normalization_3/batchnorm/mul_2:z:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization_3/batchnorm/sub?
=sequential/functional_1/batch_normalization_3/batchnorm/add_1AddV2Asequential/functional_1/batch_normalization_3/batchnorm/mul_1:z:0?sequential/functional_1/batch_normalization_3/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2?
=sequential/functional_1/batch_normalization_3/batchnorm/add_1?
3sequential/functional_1/average_pooling3d/AvgPool3D	AvgPool3DAsequential/functional_1/batch_normalization_3/batchnorm/add_1:z:0*
T0*3
_output_shapes!
:?????????*
ksize	
*
paddingVALID*
strides	
25
3sequential/functional_1/average_pooling3d/AvgPool3D?
%sequential/functional_1/flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"????   2'
%sequential/functional_1/flatten/Const?
'sequential/functional_1/flatten/ReshapeReshape<sequential/functional_1/average_pooling3d/AvgPool3D:output:0.sequential/functional_1/flatten/Const:output:0*
T0*(
_output_shapes
:??????????
2)
'sequential/functional_1/flatten/Reshape?
(sequential/functional_1/dropout/IdentityIdentity0sequential/functional_1/flatten/Reshape:output:0*
T0*(
_output_shapes
:??????????
2*
(sequential/functional_1/dropout/Identity?
4sequential/functional_1/layer1/MatMul/ReadVariableOpReadVariableOp=sequential_functional_1_layer1_matmul_readvariableop_resource* 
_output_shapes
:
?
?*
dtype026
4sequential/functional_1/layer1/MatMul/ReadVariableOp?
%sequential/functional_1/layer1/MatMulMatMul1sequential/functional_1/dropout/Identity:output:0<sequential/functional_1/layer1/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2'
%sequential/functional_1/layer1/MatMul?
5sequential/functional_1/layer1/BiasAdd/ReadVariableOpReadVariableOp>sequential_functional_1_layer1_biasadd_readvariableop_resource*
_output_shapes	
:?*
dtype027
5sequential/functional_1/layer1/BiasAdd/ReadVariableOp?
&sequential/functional_1/layer1/BiasAddBiasAdd/sequential/functional_1/layer1/MatMul:product:0=sequential/functional_1/layer1/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2(
&sequential/functional_1/layer1/BiasAdd?
"sequential/functional_1/layer1/EluElu/sequential/functional_1/layer1/BiasAdd:output:0*
T0*(
_output_shapes
:??????????2$
"sequential/functional_1/layer1/Elu?
&sequential/dense/MatMul/ReadVariableOpReadVariableOp/sequential_dense_matmul_readvariableop_resource*
_output_shapes
:	?d*
dtype02(
&sequential/dense/MatMul/ReadVariableOp?
sequential/dense/MatMulMatMul0sequential/functional_1/layer1/Elu:activations:0.sequential/dense/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????d2
sequential/dense/MatMul?
'sequential/dense/BiasAdd/ReadVariableOpReadVariableOp0sequential_dense_biasadd_readvariableop_resource*
_output_shapes
:d*
dtype02)
'sequential/dense/BiasAdd/ReadVariableOp?
sequential/dense/BiasAddBiasAdd!sequential/dense/MatMul:product:0/sequential/dense/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????d2
sequential/dense/BiasAdd?
sequential/dense/EluElu!sequential/dense/BiasAdd:output:0*
T0*'
_output_shapes
:?????????d2
sequential/dense/Elu?
(sequential/dense_1/MatMul/ReadVariableOpReadVariableOp1sequential_dense_1_matmul_readvariableop_resource*
_output_shapes

:d
*
dtype02*
(sequential/dense_1/MatMul/ReadVariableOp?
sequential/dense_1/MatMulMatMul"sequential/dense/Elu:activations:00sequential/dense_1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????
2
sequential/dense_1/MatMul?
)sequential/dense_1/BiasAdd/ReadVariableOpReadVariableOp2sequential_dense_1_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02+
)sequential/dense_1/BiasAdd/ReadVariableOp?
sequential/dense_1/BiasAddBiasAdd#sequential/dense_1/MatMul:product:01sequential/dense_1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????
2
sequential/dense_1/BiasAdd?
sequential/dense_1/EluElu#sequential/dense_1/BiasAdd:output:0*
T0*'
_output_shapes
:?????????
2
sequential/dense_1/Elu?
6sequential/functional_1/conv3d/Conv3D_1/ReadVariableOpReadVariableOp=sequential_functional_1_conv3d_conv3d_readvariableop_resource*+
_output_shapes
:?*
dtype028
6sequential/functional_1/conv3d/Conv3D_1/ReadVariableOp?
'sequential/functional_1/conv3d/Conv3D_1Conv3Dinputs_1>sequential/functional_1/conv3d/Conv3D_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2)
'sequential/functional_1/conv3d/Conv3D_1?
7sequential/functional_1/conv3d/BiasAdd_1/ReadVariableOpReadVariableOp>sequential_functional_1_conv3d_biasadd_readvariableop_resource*
_output_shapes
:*
dtype029
7sequential/functional_1/conv3d/BiasAdd_1/ReadVariableOp?
(sequential/functional_1/conv3d/BiasAdd_1BiasAdd0sequential/functional_1/conv3d/Conv3D_1:output:0?sequential/functional_1/conv3d/BiasAdd_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2*
(sequential/functional_1/conv3d/BiasAdd_1?
Fsequential/functional_1/batch_normalization/batchnorm_1/ReadVariableOpReadVariableOpMsequential_functional_1_batch_normalization_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02H
Fsequential/functional_1/batch_normalization/batchnorm_1/ReadVariableOp?
=sequential/functional_1/batch_normalization/batchnorm_1/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2?
=sequential/functional_1/batch_normalization/batchnorm_1/add/y?
;sequential/functional_1/batch_normalization/batchnorm_1/addAddV2Nsequential/functional_1/batch_normalization/batchnorm_1/ReadVariableOp:value:0Fsequential/functional_1/batch_normalization/batchnorm_1/add/y:output:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization/batchnorm_1/add?
=sequential/functional_1/batch_normalization/batchnorm_1/RsqrtRsqrt?sequential/functional_1/batch_normalization/batchnorm_1/add:z:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization/batchnorm_1/Rsqrt?
Jsequential/functional_1/batch_normalization/batchnorm_1/mul/ReadVariableOpReadVariableOpQsequential_functional_1_batch_normalization_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02L
Jsequential/functional_1/batch_normalization/batchnorm_1/mul/ReadVariableOp?
;sequential/functional_1/batch_normalization/batchnorm_1/mulMulAsequential/functional_1/batch_normalization/batchnorm_1/Rsqrt:y:0Rsequential/functional_1/batch_normalization/batchnorm_1/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization/batchnorm_1/mul?
=sequential/functional_1/batch_normalization/batchnorm_1/mul_1Mul1sequential/functional_1/conv3d/BiasAdd_1:output:0?sequential/functional_1/batch_normalization/batchnorm_1/mul:z:0*
T0*3
_output_shapes!
:?????????2?
=sequential/functional_1/batch_normalization/batchnorm_1/mul_1?
Hsequential/functional_1/batch_normalization/batchnorm_1/ReadVariableOp_1ReadVariableOpOsequential_functional_1_batch_normalization_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02J
Hsequential/functional_1/batch_normalization/batchnorm_1/ReadVariableOp_1?
=sequential/functional_1/batch_normalization/batchnorm_1/mul_2MulPsequential/functional_1/batch_normalization/batchnorm_1/ReadVariableOp_1:value:0?sequential/functional_1/batch_normalization/batchnorm_1/mul:z:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization/batchnorm_1/mul_2?
Hsequential/functional_1/batch_normalization/batchnorm_1/ReadVariableOp_2ReadVariableOpOsequential_functional_1_batch_normalization_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02J
Hsequential/functional_1/batch_normalization/batchnorm_1/ReadVariableOp_2?
;sequential/functional_1/batch_normalization/batchnorm_1/subSubPsequential/functional_1/batch_normalization/batchnorm_1/ReadVariableOp_2:value:0Asequential/functional_1/batch_normalization/batchnorm_1/mul_2:z:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization/batchnorm_1/sub?
=sequential/functional_1/batch_normalization/batchnorm_1/add_1AddV2Asequential/functional_1/batch_normalization/batchnorm_1/mul_1:z:0?sequential/functional_1/batch_normalization/batchnorm_1/sub:z:0*
T0*3
_output_shapes!
:?????????2?
=sequential/functional_1/batch_normalization/batchnorm_1/add_1?
8sequential/functional_1/conv3d_1/Conv3D_1/ReadVariableOpReadVariableOp?sequential_functional_1_conv3d_1_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02:
8sequential/functional_1/conv3d_1/Conv3D_1/ReadVariableOp?
)sequential/functional_1/conv3d_1/Conv3D_1Conv3DAsequential/functional_1/batch_normalization/batchnorm_1/add_1:z:0@sequential/functional_1/conv3d_1/Conv3D_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2+
)sequential/functional_1/conv3d_1/Conv3D_1?
9sequential/functional_1/conv3d_1/BiasAdd_1/ReadVariableOpReadVariableOp@sequential_functional_1_conv3d_1_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02;
9sequential/functional_1/conv3d_1/BiasAdd_1/ReadVariableOp?
*sequential/functional_1/conv3d_1/BiasAdd_1BiasAdd2sequential/functional_1/conv3d_1/Conv3D_1:output:0Asequential/functional_1/conv3d_1/BiasAdd_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2,
*sequential/functional_1/conv3d_1/BiasAdd_1?
&sequential/functional_1/conv3d_1/Elu_1Elu3sequential/functional_1/conv3d_1/BiasAdd_1:output:0*
T0*3
_output_shapes!
:?????????2(
&sequential/functional_1/conv3d_1/Elu_1?
Hsequential/functional_1/batch_normalization_1/batchnorm_1/ReadVariableOpReadVariableOpOsequential_functional_1_batch_normalization_1_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02J
Hsequential/functional_1/batch_normalization_1/batchnorm_1/ReadVariableOp?
?sequential/functional_1/batch_normalization_1/batchnorm_1/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2A
?sequential/functional_1/batch_normalization_1/batchnorm_1/add/y?
=sequential/functional_1/batch_normalization_1/batchnorm_1/addAddV2Psequential/functional_1/batch_normalization_1/batchnorm_1/ReadVariableOp:value:0Hsequential/functional_1/batch_normalization_1/batchnorm_1/add/y:output:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_1/batchnorm_1/add?
?sequential/functional_1/batch_normalization_1/batchnorm_1/RsqrtRsqrtAsequential/functional_1/batch_normalization_1/batchnorm_1/add:z:0*
T0*
_output_shapes
:2A
?sequential/functional_1/batch_normalization_1/batchnorm_1/Rsqrt?
Lsequential/functional_1/batch_normalization_1/batchnorm_1/mul/ReadVariableOpReadVariableOpSsequential_functional_1_batch_normalization_1_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02N
Lsequential/functional_1/batch_normalization_1/batchnorm_1/mul/ReadVariableOp?
=sequential/functional_1/batch_normalization_1/batchnorm_1/mulMulCsequential/functional_1/batch_normalization_1/batchnorm_1/Rsqrt:y:0Tsequential/functional_1/batch_normalization_1/batchnorm_1/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_1/batchnorm_1/mul?
?sequential/functional_1/batch_normalization_1/batchnorm_1/mul_1Mul4sequential/functional_1/conv3d_1/Elu_1:activations:0Asequential/functional_1/batch_normalization_1/batchnorm_1/mul:z:0*
T0*3
_output_shapes!
:?????????2A
?sequential/functional_1/batch_normalization_1/batchnorm_1/mul_1?
Jsequential/functional_1/batch_normalization_1/batchnorm_1/ReadVariableOp_1ReadVariableOpQsequential_functional_1_batch_normalization_1_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02L
Jsequential/functional_1/batch_normalization_1/batchnorm_1/ReadVariableOp_1?
?sequential/functional_1/batch_normalization_1/batchnorm_1/mul_2MulRsequential/functional_1/batch_normalization_1/batchnorm_1/ReadVariableOp_1:value:0Asequential/functional_1/batch_normalization_1/batchnorm_1/mul:z:0*
T0*
_output_shapes
:2A
?sequential/functional_1/batch_normalization_1/batchnorm_1/mul_2?
Jsequential/functional_1/batch_normalization_1/batchnorm_1/ReadVariableOp_2ReadVariableOpQsequential_functional_1_batch_normalization_1_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02L
Jsequential/functional_1/batch_normalization_1/batchnorm_1/ReadVariableOp_2?
=sequential/functional_1/batch_normalization_1/batchnorm_1/subSubRsequential/functional_1/batch_normalization_1/batchnorm_1/ReadVariableOp_2:value:0Csequential/functional_1/batch_normalization_1/batchnorm_1/mul_2:z:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_1/batchnorm_1/sub?
?sequential/functional_1/batch_normalization_1/batchnorm_1/add_1AddV2Csequential/functional_1/batch_normalization_1/batchnorm_1/mul_1:z:0Asequential/functional_1/batch_normalization_1/batchnorm_1/sub:z:0*
T0*3
_output_shapes!
:?????????2A
?sequential/functional_1/batch_normalization_1/batchnorm_1/add_1?
8sequential/functional_1/conv3d_2/Conv3D_1/ReadVariableOpReadVariableOp?sequential_functional_1_conv3d_2_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02:
8sequential/functional_1/conv3d_2/Conv3D_1/ReadVariableOp?
)sequential/functional_1/conv3d_2/Conv3D_1Conv3DCsequential/functional_1/batch_normalization_1/batchnorm_1/add_1:z:0@sequential/functional_1/conv3d_2/Conv3D_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2+
)sequential/functional_1/conv3d_2/Conv3D_1?
9sequential/functional_1/conv3d_2/BiasAdd_1/ReadVariableOpReadVariableOp@sequential_functional_1_conv3d_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02;
9sequential/functional_1/conv3d_2/BiasAdd_1/ReadVariableOp?
*sequential/functional_1/conv3d_2/BiasAdd_1BiasAdd2sequential/functional_1/conv3d_2/Conv3D_1:output:0Asequential/functional_1/conv3d_2/BiasAdd_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2,
*sequential/functional_1/conv3d_2/BiasAdd_1?
&sequential/functional_1/conv3d_2/Elu_1Elu3sequential/functional_1/conv3d_2/BiasAdd_1:output:0*
T0*3
_output_shapes!
:?????????2(
&sequential/functional_1/conv3d_2/Elu_1?
Hsequential/functional_1/batch_normalization_2/batchnorm_1/ReadVariableOpReadVariableOpOsequential_functional_1_batch_normalization_2_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02J
Hsequential/functional_1/batch_normalization_2/batchnorm_1/ReadVariableOp?
?sequential/functional_1/batch_normalization_2/batchnorm_1/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2A
?sequential/functional_1/batch_normalization_2/batchnorm_1/add/y?
=sequential/functional_1/batch_normalization_2/batchnorm_1/addAddV2Psequential/functional_1/batch_normalization_2/batchnorm_1/ReadVariableOp:value:0Hsequential/functional_1/batch_normalization_2/batchnorm_1/add/y:output:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_2/batchnorm_1/add?
?sequential/functional_1/batch_normalization_2/batchnorm_1/RsqrtRsqrtAsequential/functional_1/batch_normalization_2/batchnorm_1/add:z:0*
T0*
_output_shapes
:2A
?sequential/functional_1/batch_normalization_2/batchnorm_1/Rsqrt?
Lsequential/functional_1/batch_normalization_2/batchnorm_1/mul/ReadVariableOpReadVariableOpSsequential_functional_1_batch_normalization_2_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02N
Lsequential/functional_1/batch_normalization_2/batchnorm_1/mul/ReadVariableOp?
=sequential/functional_1/batch_normalization_2/batchnorm_1/mulMulCsequential/functional_1/batch_normalization_2/batchnorm_1/Rsqrt:y:0Tsequential/functional_1/batch_normalization_2/batchnorm_1/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_2/batchnorm_1/mul?
?sequential/functional_1/batch_normalization_2/batchnorm_1/mul_1Mul4sequential/functional_1/conv3d_2/Elu_1:activations:0Asequential/functional_1/batch_normalization_2/batchnorm_1/mul:z:0*
T0*3
_output_shapes!
:?????????2A
?sequential/functional_1/batch_normalization_2/batchnorm_1/mul_1?
Jsequential/functional_1/batch_normalization_2/batchnorm_1/ReadVariableOp_1ReadVariableOpQsequential_functional_1_batch_normalization_2_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02L
Jsequential/functional_1/batch_normalization_2/batchnorm_1/ReadVariableOp_1?
?sequential/functional_1/batch_normalization_2/batchnorm_1/mul_2MulRsequential/functional_1/batch_normalization_2/batchnorm_1/ReadVariableOp_1:value:0Asequential/functional_1/batch_normalization_2/batchnorm_1/mul:z:0*
T0*
_output_shapes
:2A
?sequential/functional_1/batch_normalization_2/batchnorm_1/mul_2?
Jsequential/functional_1/batch_normalization_2/batchnorm_1/ReadVariableOp_2ReadVariableOpQsequential_functional_1_batch_normalization_2_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02L
Jsequential/functional_1/batch_normalization_2/batchnorm_1/ReadVariableOp_2?
=sequential/functional_1/batch_normalization_2/batchnorm_1/subSubRsequential/functional_1/batch_normalization_2/batchnorm_1/ReadVariableOp_2:value:0Csequential/functional_1/batch_normalization_2/batchnorm_1/mul_2:z:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_2/batchnorm_1/sub?
?sequential/functional_1/batch_normalization_2/batchnorm_1/add_1AddV2Csequential/functional_1/batch_normalization_2/batchnorm_1/mul_1:z:0Asequential/functional_1/batch_normalization_2/batchnorm_1/sub:z:0*
T0*3
_output_shapes!
:?????????2A
?sequential/functional_1/batch_normalization_2/batchnorm_1/add_1?
8sequential/functional_1/conv3d_3/Conv3D_1/ReadVariableOpReadVariableOp?sequential_functional_1_conv3d_3_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02:
8sequential/functional_1/conv3d_3/Conv3D_1/ReadVariableOp?
)sequential/functional_1/conv3d_3/Conv3D_1Conv3DCsequential/functional_1/batch_normalization_2/batchnorm_1/add_1:z:0@sequential/functional_1/conv3d_3/Conv3D_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2+
)sequential/functional_1/conv3d_3/Conv3D_1?
9sequential/functional_1/conv3d_3/BiasAdd_1/ReadVariableOpReadVariableOp@sequential_functional_1_conv3d_3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02;
9sequential/functional_1/conv3d_3/BiasAdd_1/ReadVariableOp?
*sequential/functional_1/conv3d_3/BiasAdd_1BiasAdd2sequential/functional_1/conv3d_3/Conv3D_1:output:0Asequential/functional_1/conv3d_3/BiasAdd_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2,
*sequential/functional_1/conv3d_3/BiasAdd_1?
&sequential/functional_1/conv3d_3/Elu_1Elu3sequential/functional_1/conv3d_3/BiasAdd_1:output:0*
T0*3
_output_shapes!
:?????????2(
&sequential/functional_1/conv3d_3/Elu_1?
Hsequential/functional_1/batch_normalization_3/batchnorm_1/ReadVariableOpReadVariableOpOsequential_functional_1_batch_normalization_3_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02J
Hsequential/functional_1/batch_normalization_3/batchnorm_1/ReadVariableOp?
?sequential/functional_1/batch_normalization_3/batchnorm_1/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2A
?sequential/functional_1/batch_normalization_3/batchnorm_1/add/y?
=sequential/functional_1/batch_normalization_3/batchnorm_1/addAddV2Psequential/functional_1/batch_normalization_3/batchnorm_1/ReadVariableOp:value:0Hsequential/functional_1/batch_normalization_3/batchnorm_1/add/y:output:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_3/batchnorm_1/add?
?sequential/functional_1/batch_normalization_3/batchnorm_1/RsqrtRsqrtAsequential/functional_1/batch_normalization_3/batchnorm_1/add:z:0*
T0*
_output_shapes
:2A
?sequential/functional_1/batch_normalization_3/batchnorm_1/Rsqrt?
Lsequential/functional_1/batch_normalization_3/batchnorm_1/mul/ReadVariableOpReadVariableOpSsequential_functional_1_batch_normalization_3_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02N
Lsequential/functional_1/batch_normalization_3/batchnorm_1/mul/ReadVariableOp?
=sequential/functional_1/batch_normalization_3/batchnorm_1/mulMulCsequential/functional_1/batch_normalization_3/batchnorm_1/Rsqrt:y:0Tsequential/functional_1/batch_normalization_3/batchnorm_1/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_3/batchnorm_1/mul?
?sequential/functional_1/batch_normalization_3/batchnorm_1/mul_1Mul4sequential/functional_1/conv3d_3/Elu_1:activations:0Asequential/functional_1/batch_normalization_3/batchnorm_1/mul:z:0*
T0*3
_output_shapes!
:?????????2A
?sequential/functional_1/batch_normalization_3/batchnorm_1/mul_1?
Jsequential/functional_1/batch_normalization_3/batchnorm_1/ReadVariableOp_1ReadVariableOpQsequential_functional_1_batch_normalization_3_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02L
Jsequential/functional_1/batch_normalization_3/batchnorm_1/ReadVariableOp_1?
?sequential/functional_1/batch_normalization_3/batchnorm_1/mul_2MulRsequential/functional_1/batch_normalization_3/batchnorm_1/ReadVariableOp_1:value:0Asequential/functional_1/batch_normalization_3/batchnorm_1/mul:z:0*
T0*
_output_shapes
:2A
?sequential/functional_1/batch_normalization_3/batchnorm_1/mul_2?
Jsequential/functional_1/batch_normalization_3/batchnorm_1/ReadVariableOp_2ReadVariableOpQsequential_functional_1_batch_normalization_3_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02L
Jsequential/functional_1/batch_normalization_3/batchnorm_1/ReadVariableOp_2?
=sequential/functional_1/batch_normalization_3/batchnorm_1/subSubRsequential/functional_1/batch_normalization_3/batchnorm_1/ReadVariableOp_2:value:0Csequential/functional_1/batch_normalization_3/batchnorm_1/mul_2:z:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_3/batchnorm_1/sub?
?sequential/functional_1/batch_normalization_3/batchnorm_1/add_1AddV2Csequential/functional_1/batch_normalization_3/batchnorm_1/mul_1:z:0Asequential/functional_1/batch_normalization_3/batchnorm_1/sub:z:0*
T0*3
_output_shapes!
:?????????2A
?sequential/functional_1/batch_normalization_3/batchnorm_1/add_1?
5sequential/functional_1/average_pooling3d/AvgPool3D_1	AvgPool3DCsequential/functional_1/batch_normalization_3/batchnorm_1/add_1:z:0*
T0*3
_output_shapes!
:?????????*
ksize	
*
paddingVALID*
strides	
27
5sequential/functional_1/average_pooling3d/AvgPool3D_1?
'sequential/functional_1/flatten/Const_1Const*
_output_shapes
:*
dtype0*
valueB"????   2)
'sequential/functional_1/flatten/Const_1?
)sequential/functional_1/flatten/Reshape_1Reshape>sequential/functional_1/average_pooling3d/AvgPool3D_1:output:00sequential/functional_1/flatten/Const_1:output:0*
T0*(
_output_shapes
:??????????
2+
)sequential/functional_1/flatten/Reshape_1?
*sequential/functional_1/dropout/Identity_1Identity2sequential/functional_1/flatten/Reshape_1:output:0*
T0*(
_output_shapes
:??????????
2,
*sequential/functional_1/dropout/Identity_1?
6sequential/functional_1/layer1/MatMul_1/ReadVariableOpReadVariableOp=sequential_functional_1_layer1_matmul_readvariableop_resource* 
_output_shapes
:
?
?*
dtype028
6sequential/functional_1/layer1/MatMul_1/ReadVariableOp?
'sequential/functional_1/layer1/MatMul_1MatMul3sequential/functional_1/dropout/Identity_1:output:0>sequential/functional_1/layer1/MatMul_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2)
'sequential/functional_1/layer1/MatMul_1?
7sequential/functional_1/layer1/BiasAdd_1/ReadVariableOpReadVariableOp>sequential_functional_1_layer1_biasadd_readvariableop_resource*
_output_shapes	
:?*
dtype029
7sequential/functional_1/layer1/BiasAdd_1/ReadVariableOp?
(sequential/functional_1/layer1/BiasAdd_1BiasAdd1sequential/functional_1/layer1/MatMul_1:product:0?sequential/functional_1/layer1/BiasAdd_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2*
(sequential/functional_1/layer1/BiasAdd_1?
$sequential/functional_1/layer1/Elu_1Elu1sequential/functional_1/layer1/BiasAdd_1:output:0*
T0*(
_output_shapes
:??????????2&
$sequential/functional_1/layer1/Elu_1?
(sequential/dense/MatMul_1/ReadVariableOpReadVariableOp/sequential_dense_matmul_readvariableop_resource*
_output_shapes
:	?d*
dtype02*
(sequential/dense/MatMul_1/ReadVariableOp?
sequential/dense/MatMul_1MatMul2sequential/functional_1/layer1/Elu_1:activations:00sequential/dense/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????d2
sequential/dense/MatMul_1?
)sequential/dense/BiasAdd_1/ReadVariableOpReadVariableOp0sequential_dense_biasadd_readvariableop_resource*
_output_shapes
:d*
dtype02+
)sequential/dense/BiasAdd_1/ReadVariableOp?
sequential/dense/BiasAdd_1BiasAdd#sequential/dense/MatMul_1:product:01sequential/dense/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????d2
sequential/dense/BiasAdd_1?
sequential/dense/Elu_1Elu#sequential/dense/BiasAdd_1:output:0*
T0*'
_output_shapes
:?????????d2
sequential/dense/Elu_1?
*sequential/dense_1/MatMul_1/ReadVariableOpReadVariableOp1sequential_dense_1_matmul_readvariableop_resource*
_output_shapes

:d
*
dtype02,
*sequential/dense_1/MatMul_1/ReadVariableOp?
sequential/dense_1/MatMul_1MatMul$sequential/dense/Elu_1:activations:02sequential/dense_1/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????
2
sequential/dense_1/MatMul_1?
+sequential/dense_1/BiasAdd_1/ReadVariableOpReadVariableOp2sequential_dense_1_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02-
+sequential/dense_1/BiasAdd_1/ReadVariableOp?
sequential/dense_1/BiasAdd_1BiasAdd%sequential/dense_1/MatMul_1:product:03sequential/dense_1/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????
2
sequential/dense_1/BiasAdd_1?
sequential/dense_1/Elu_1Elu%sequential/dense_1/BiasAdd_1:output:0*
T0*'
_output_shapes
:?????????
2
sequential/dense_1/Elu_1?

lambda/subSub$sequential/dense_1/Elu:activations:0&sequential/dense_1/Elu_1:activations:0*
T0*'
_output_shapes
:?????????
2

lambda/suba

lambda/AbsAbslambda/sub:z:0*
T0*'
_output_shapes
:?????????
2

lambda/Abst
concatenate/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
concatenate/concat/axis?
concatenate/concatConcatV2lambda/Abs:y:0inputs_2 concatenate/concat/axis:output:0*
N*
T0*'
_output_shapes
:?????????2
concatenate/concat?
dense_2/MatMul/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes

:*
dtype02
dense_2/MatMul/ReadVariableOp?
dense_2/MatMulMatMulconcatenate/concat:output:0%dense_2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
dense_2/MatMul?
dense_2/BiasAdd/ReadVariableOpReadVariableOp'dense_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02 
dense_2/BiasAdd/ReadVariableOp?
dense_2/BiasAddBiasAdddense_2/MatMul:product:0&dense_2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
dense_2/BiasAdd?
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOp/sequential_dense_matmul_readvariableop_resource*
_output_shapes
:	?d*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOp?
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?d2!
dense/kernel/Regularizer/Square?
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const?
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/Sum?
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2 
dense/kernel/Regularizer/mul/x?
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/mul?
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp1sequential_dense_1_matmul_readvariableop_resource*
_output_shapes

:d
*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp?
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d
2#
!dense_1/kernel/Regularizer/Square?
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const?
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/Sum?
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2"
 dense_1/kernel/Regularizer/mul/x?
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/mull
IdentityIdentitydense_2/BiasAdd:output:0*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????:??????????:?????????
:::::::::::::::::::::::::::::::::^ Z
4
_output_shapes"
 :??????????
"
_user_specified_name
inputs/0:^Z
4
_output_shapes"
 :??????????
"
_user_specified_name
inputs/1:QM
'
_output_shapes
:?????????

"
_user_specified_name
inputs/2
?
?
7__inference_batch_normalization_3_layer_call_fn_1227645

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_12232032
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????::::22
StatefulPartitionedCallStatefulPartitionedCall:[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?
}
(__inference_layer1_layer_call_fn_1227703

inputs
unknown
	unknown_0
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *L
fGRE
C__inference_layer1_layer_call_and_return_conditional_losses_12232952
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*/
_input_shapes
:??????????
::22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:??????????

 
_user_specified_nameinputs
?
?
7__inference_batch_normalization_3_layer_call_fn_1227563

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *N
_output_shapes<
::8????????????????????????????????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_12227412
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????::::22
StatefulPartitionedCallStatefulPartitionedCall:v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
?
?
C__inference_conv3d_layer_call_and_return_conditional_losses_1226920

inputs"
conv3d_readvariableop_resource#
biasadd_readvariableop_resource
identity??
Conv3D/ReadVariableOpReadVariableOpconv3d_readvariableop_resource*+
_output_shapes
:?*
dtype02
Conv3D/ReadVariableOp?
Conv3DConv3DinputsConv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
Conv3D?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddConv3D:output:0BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2	
BiasAddp
IdentityIdentityBiasAdd:output:0*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*;
_input_shapes*
(:??????????:::\ X
4
_output_shapes"
 :??????????
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_1223203

inputs%
!batchnorm_readvariableop_resource)
%batchnorm_mul_readvariableop_resource'
#batchnorm_readvariableop_1_resource'
#batchnorm_readvariableop_2_resource
identity??
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/mul_1?
batchnorm/ReadVariableOp_1ReadVariableOp#batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_1?
batchnorm/mul_2Mul"batchnorm/ReadVariableOp_1:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOp_2ReadVariableOp#batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_2?
batchnorm/subSub"batchnorm/ReadVariableOp_2:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/add_1s
IdentityIdentitybatchnorm/add_1:z:0*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????:::::[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?
t
H__inference_concatenate_layer_call_and_return_conditional_losses_1226402
inputs_0
inputs_1
identity\
concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
concat/axis?
concatConcatV2inputs_0inputs_1concat/axis:output:0*
N*
T0*'
_output_shapes
:?????????2
concatc
IdentityIdentityconcat:output:0*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*9
_input_shapes(
&:?????????
:?????????
:Q M
'
_output_shapes
:?????????

"
_user_specified_name
inputs/0:QM
'
_output_shapes
:?????????

"
_user_specified_name
inputs/1
?B
?	
I__inference_functional_1_layer_call_and_return_conditional_losses_1223380
input_1
conv3d_1223315
conv3d_1223317
batch_normalization_1223320
batch_normalization_1223322
batch_normalization_1223324
batch_normalization_1223326
conv3d_1_1223329
conv3d_1_1223331!
batch_normalization_1_1223334!
batch_normalization_1_1223336!
batch_normalization_1_1223338!
batch_normalization_1_1223340
conv3d_2_1223343
conv3d_2_1223345!
batch_normalization_2_1223348!
batch_normalization_2_1223350!
batch_normalization_2_1223352!
batch_normalization_2_1223354
conv3d_3_1223357
conv3d_3_1223359!
batch_normalization_3_1223362!
batch_normalization_3_1223364!
batch_normalization_3_1223366!
batch_normalization_3_1223368
layer1_1223374
layer1_1223376
identity??+batch_normalization/StatefulPartitionedCall?-batch_normalization_1/StatefulPartitionedCall?-batch_normalization_2/StatefulPartitionedCall?-batch_normalization_3/StatefulPartitionedCall?conv3d/StatefulPartitionedCall? conv3d_1/StatefulPartitionedCall? conv3d_2/StatefulPartitionedCall? conv3d_3/StatefulPartitionedCall?layer1/StatefulPartitionedCall?
conv3d/StatefulPartitionedCallStatefulPartitionedCallinput_1conv3d_1223315conv3d_1223317*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *L
fGRE
C__inference_conv3d_layer_call_and_return_conditional_losses_12227782 
conv3d/StatefulPartitionedCall?
+batch_normalization/StatefulPartitionedCallStatefulPartitionedCall'conv3d/StatefulPartitionedCall:output:0batch_normalization_1223320batch_normalization_1223322batch_normalization_1223324batch_normalization_1223326*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *Y
fTRR
P__inference_batch_normalization_layer_call_and_return_conditional_losses_12228492-
+batch_normalization/StatefulPartitionedCall?
 conv3d_1/StatefulPartitionedCallStatefulPartitionedCall4batch_normalization/StatefulPartitionedCall:output:0conv3d_1_1223329conv3d_1_1223331*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv3d_1_layer_call_and_return_conditional_losses_12228962"
 conv3d_1/StatefulPartitionedCall?
-batch_normalization_1/StatefulPartitionedCallStatefulPartitionedCall)conv3d_1/StatefulPartitionedCall:output:0batch_normalization_1_1223334batch_normalization_1_1223336batch_normalization_1_1223338batch_normalization_1_1223340*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_12229672/
-batch_normalization_1/StatefulPartitionedCall?
 conv3d_2/StatefulPartitionedCallStatefulPartitionedCall6batch_normalization_1/StatefulPartitionedCall:output:0conv3d_2_1223343conv3d_2_1223345*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv3d_2_layer_call_and_return_conditional_losses_12230142"
 conv3d_2/StatefulPartitionedCall?
-batch_normalization_2/StatefulPartitionedCallStatefulPartitionedCall)conv3d_2/StatefulPartitionedCall:output:0batch_normalization_2_1223348batch_normalization_2_1223350batch_normalization_2_1223352batch_normalization_2_1223354*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_12230852/
-batch_normalization_2/StatefulPartitionedCall?
 conv3d_3/StatefulPartitionedCallStatefulPartitionedCall6batch_normalization_2/StatefulPartitionedCall:output:0conv3d_3_1223357conv3d_3_1223359*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv3d_3_layer_call_and_return_conditional_losses_12231322"
 conv3d_3/StatefulPartitionedCall?
-batch_normalization_3/StatefulPartitionedCallStatefulPartitionedCall)conv3d_3/StatefulPartitionedCall:output:0batch_normalization_3_1223362batch_normalization_3_1223364batch_normalization_3_1223366batch_normalization_3_1223368*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_12232032/
-batch_normalization_3/StatefulPartitionedCall?
!average_pooling3d/PartitionedCallPartitionedCall6batch_normalization_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *W
fRRP
N__inference_average_pooling3d_layer_call_and_return_conditional_losses_12227582#
!average_pooling3d/PartitionedCall?
flatten/PartitionedCallPartitionedCall*average_pooling3d/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????
* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_flatten_layer_call_and_return_conditional_losses_12232462
flatten/PartitionedCall?
dropout/PartitionedCallPartitionedCall flatten/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????
* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dropout_layer_call_and_return_conditional_losses_12232712
dropout/PartitionedCall?
layer1/StatefulPartitionedCallStatefulPartitionedCall dropout/PartitionedCall:output:0layer1_1223374layer1_1223376*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *L
fGRE
C__inference_layer1_layer_call_and_return_conditional_losses_12232952 
layer1/StatefulPartitionedCall?
IdentityIdentity'layer1/StatefulPartitionedCall:output:0,^batch_normalization/StatefulPartitionedCall.^batch_normalization_1/StatefulPartitionedCall.^batch_normalization_2/StatefulPartitionedCall.^batch_normalization_3/StatefulPartitionedCall^conv3d/StatefulPartitionedCall!^conv3d_1/StatefulPartitionedCall!^conv3d_2/StatefulPartitionedCall!^conv3d_3/StatefulPartitionedCall^layer1/StatefulPartitionedCall*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????::::::::::::::::::::::::::2Z
+batch_normalization/StatefulPartitionedCall+batch_normalization/StatefulPartitionedCall2^
-batch_normalization_1/StatefulPartitionedCall-batch_normalization_1/StatefulPartitionedCall2^
-batch_normalization_2/StatefulPartitionedCall-batch_normalization_2/StatefulPartitionedCall2^
-batch_normalization_3/StatefulPartitionedCall-batch_normalization_3/StatefulPartitionedCall2@
conv3d/StatefulPartitionedCallconv3d/StatefulPartitionedCall2D
 conv3d_1/StatefulPartitionedCall conv3d_1/StatefulPartitionedCall2D
 conv3d_2/StatefulPartitionedCall conv3d_2/StatefulPartitionedCall2D
 conv3d_3/StatefulPartitionedCall conv3d_3/StatefulPartitionedCall2@
layer1/StatefulPartitionedCalllayer1/StatefulPartitionedCall:] Y
4
_output_shapes"
 :??????????
!
_user_specified_name	input_1
?
?
.__inference_functional_1_layer_call_fn_1223506
input_1
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16

unknown_17

unknown_18

unknown_19

unknown_20

unknown_21

unknown_22

unknown_23

unknown_24
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19
unknown_20
unknown_21
unknown_22
unknown_23
unknown_24*&
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*4
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *R
fMRK
I__inference_functional_1_layer_call_and_return_conditional_losses_12234512
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????::::::::::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:] Y
4
_output_shapes"
 :??????????
!
_user_specified_name	input_1
?+
?
G__inference_sequential_layer_call_and_return_conditional_losses_1224042

inputs
functional_1_1223966
functional_1_1223968
functional_1_1223970
functional_1_1223972
functional_1_1223974
functional_1_1223976
functional_1_1223978
functional_1_1223980
functional_1_1223982
functional_1_1223984
functional_1_1223986
functional_1_1223988
functional_1_1223990
functional_1_1223992
functional_1_1223994
functional_1_1223996
functional_1_1223998
functional_1_1224000
functional_1_1224002
functional_1_1224004
functional_1_1224006
functional_1_1224008
functional_1_1224010
functional_1_1224012
functional_1_1224014
functional_1_1224016
dense_1224019
dense_1224021
dense_1_1224024
dense_1_1224026
identity??dense/StatefulPartitionedCall?dense_1/StatefulPartitionedCall?$functional_1/StatefulPartitionedCall?
$functional_1/StatefulPartitionedCallStatefulPartitionedCallinputsfunctional_1_1223966functional_1_1223968functional_1_1223970functional_1_1223972functional_1_1223974functional_1_1223976functional_1_1223978functional_1_1223980functional_1_1223982functional_1_1223984functional_1_1223986functional_1_1223988functional_1_1223990functional_1_1223992functional_1_1223994functional_1_1223996functional_1_1223998functional_1_1224000functional_1_1224002functional_1_1224004functional_1_1224006functional_1_1224008functional_1_1224010functional_1_1224012functional_1_1224014functional_1_1224016*&
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*4
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *R
fMRK
I__inference_functional_1_layer_call_and_return_conditional_losses_12234512&
$functional_1/StatefulPartitionedCall?
dense/StatefulPartitionedCallStatefulPartitionedCall-functional_1/StatefulPartitionedCall:output:0dense_1224019dense_1224021*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????d*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *K
fFRD
B__inference_dense_layer_call_and_return_conditional_losses_12238192
dense/StatefulPartitionedCall?
dense_1/StatefulPartitionedCallStatefulPartitionedCall&dense/StatefulPartitionedCall:output:0dense_1_1224024dense_1_1224026*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_1_layer_call_and_return_conditional_losses_12238522!
dense_1/StatefulPartitionedCall?
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_1224019*
_output_shapes
:	?d*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOp?
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?d2!
dense/kernel/Regularizer/Square?
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const?
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/Sum?
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2 
dense/kernel/Regularizer/mul/x?
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/mul?
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_1224024*
_output_shapes

:d
*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp?
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d
2#
!dense_1/kernel/Regularizer/Square?
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const?
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/Sum?
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2"
 dense_1/kernel/Regularizer/mul/x?
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/mul?
IdentityIdentity(dense_1/StatefulPartitionedCall:output:0^dense/StatefulPartitionedCall ^dense_1/StatefulPartitionedCall%^functional_1/StatefulPartitionedCall*
T0*'
_output_shapes
:?????????
2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????::::::::::::::::::::::::::::::2>
dense/StatefulPartitionedCalldense/StatefulPartitionedCall2B
dense_1/StatefulPartitionedCalldense_1/StatefulPartitionedCall2L
$functional_1/StatefulPartitionedCall$functional_1/StatefulPartitionedCall:\ X
4
_output_shapes"
 :??????????
 
_user_specified_nameinputs
?
b
D__inference_dropout_layer_call_and_return_conditional_losses_1227673

inputs

identity_1[
IdentityIdentityinputs*
T0*(
_output_shapes
:??????????
2

Identityj

Identity_1IdentityIdentity:output:0*
T0*(
_output_shapes
:??????????
2

Identity_1"!

identity_1Identity_1:output:0*'
_input_shapes
:??????????
:P L
(
_output_shapes
:??????????

 
_user_specified_nameinputs
?9
?
I__inference_functional_3_layer_call_and_return_conditional_losses_1224561
input_1
input_2
input_3
sequential_1224385
sequential_1224387
sequential_1224389
sequential_1224391
sequential_1224393
sequential_1224395
sequential_1224397
sequential_1224399
sequential_1224401
sequential_1224403
sequential_1224405
sequential_1224407
sequential_1224409
sequential_1224411
sequential_1224413
sequential_1224415
sequential_1224417
sequential_1224419
sequential_1224421
sequential_1224423
sequential_1224425
sequential_1224427
sequential_1224429
sequential_1224431
sequential_1224433
sequential_1224435
sequential_1224437
sequential_1224439
sequential_1224441
sequential_1224443
dense_2_1224543
dense_2_1224545
identity??dense_2/StatefulPartitionedCall?"sequential/StatefulPartitionedCall?$sequential/StatefulPartitionedCall_1?
"sequential/StatefulPartitionedCallStatefulPartitionedCallinput_1sequential_1224385sequential_1224387sequential_1224389sequential_1224391sequential_1224393sequential_1224395sequential_1224397sequential_1224399sequential_1224401sequential_1224403sequential_1224405sequential_1224407sequential_1224409sequential_1224411sequential_1224413sequential_1224415sequential_1224417sequential_1224419sequential_1224421sequential_1224423sequential_1224425sequential_1224427sequential_1224429sequential_1224431sequential_1224433sequential_1224435sequential_1224437sequential_1224439sequential_1224441sequential_1224443**
Tin#
!2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
*8
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *P
fKRI
G__inference_sequential_layer_call_and_return_conditional_losses_12240422$
"sequential/StatefulPartitionedCall?
$sequential/StatefulPartitionedCall_1StatefulPartitionedCallinput_2sequential_1224385sequential_1224387sequential_1224389sequential_1224391sequential_1224393sequential_1224395sequential_1224397sequential_1224399sequential_1224401sequential_1224403sequential_1224405sequential_1224407sequential_1224409sequential_1224411sequential_1224413sequential_1224415sequential_1224417sequential_1224419sequential_1224421sequential_1224423sequential_1224425sequential_1224427sequential_1224429sequential_1224431sequential_1224433sequential_1224435sequential_1224437sequential_1224439sequential_1224441sequential_1224443#^sequential/StatefulPartitionedCall**
Tin#
!2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
*8
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *P
fKRI
G__inference_sequential_layer_call_and_return_conditional_losses_12240422&
$sequential/StatefulPartitionedCall_1?
lambda/PartitionedCallPartitionedCall+sequential/StatefulPartitionedCall:output:0-sequential/StatefulPartitionedCall_1:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *L
fGRE
C__inference_lambda_layer_call_and_return_conditional_losses_12244842
lambda/PartitionedCall?
concatenate/PartitionedCallPartitionedCalllambda/PartitionedCall:output:0input_3*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *Q
fLRJ
H__inference_concatenate_layer_call_and_return_conditional_losses_12245132
concatenate/PartitionedCall?
dense_2/StatefulPartitionedCallStatefulPartitionedCall$concatenate/PartitionedCall:output:0dense_2_1224543dense_2_1224545*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_2_layer_call_and_return_conditional_losses_12245322!
dense_2/StatefulPartitionedCall?
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_1224437*
_output_shapes
:	?d*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOp?
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?d2!
dense/kernel/Regularizer/Square?
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const?
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/Sum?
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2 
dense/kernel/Regularizer/mul/x?
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/mul?
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_1224441*
_output_shapes

:d
*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp?
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d
2#
!dense_1/kernel/Regularizer/Square?
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const?
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/Sum?
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2"
 dense_1/kernel/Regularizer/mul/x?
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/mul?
IdentityIdentity(dense_2/StatefulPartitionedCall:output:0 ^dense_2/StatefulPartitionedCall#^sequential/StatefulPartitionedCall%^sequential/StatefulPartitionedCall_1*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????:??????????:?????????
::::::::::::::::::::::::::::::::2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2H
"sequential/StatefulPartitionedCall"sequential/StatefulPartitionedCall2L
$sequential/StatefulPartitionedCall_1$sequential/StatefulPartitionedCall_1:] Y
4
_output_shapes"
 :??????????
!
_user_specified_name	input_1:]Y
4
_output_shapes"
 :??????????
!
_user_specified_name	input_2:PL
'
_output_shapes
:?????????

!
_user_specified_name	input_3
?
?
7__inference_batch_normalization_2_layer_call_fn_1227448

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_12230652
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????::::22
StatefulPartitionedCallStatefulPartitionedCall:[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?
?
C__inference_layer1_layer_call_and_return_conditional_losses_1227694

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity??
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
?
?*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2
MatMul?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:?*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2	
BiasAddV
EluEluBiasAdd:output:0*
T0*(
_output_shapes
:??????????2
Eluf
IdentityIdentityElu:activations:0*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*/
_input_shapes
:??????????
:::P L
(
_output_shapes
:??????????

 
_user_specified_nameinputs
?9
?
I__inference_functional_3_layer_call_and_return_conditional_losses_1224988

inputs
inputs_1
inputs_2
sequential_1224876
sequential_1224878
sequential_1224880
sequential_1224882
sequential_1224884
sequential_1224886
sequential_1224888
sequential_1224890
sequential_1224892
sequential_1224894
sequential_1224896
sequential_1224898
sequential_1224900
sequential_1224902
sequential_1224904
sequential_1224906
sequential_1224908
sequential_1224910
sequential_1224912
sequential_1224914
sequential_1224916
sequential_1224918
sequential_1224920
sequential_1224922
sequential_1224924
sequential_1224926
sequential_1224928
sequential_1224930
sequential_1224932
sequential_1224934
dense_2_1224970
dense_2_1224972
identity??dense_2/StatefulPartitionedCall?"sequential/StatefulPartitionedCall?$sequential/StatefulPartitionedCall_1?
"sequential/StatefulPartitionedCallStatefulPartitionedCallinputssequential_1224876sequential_1224878sequential_1224880sequential_1224882sequential_1224884sequential_1224886sequential_1224888sequential_1224890sequential_1224892sequential_1224894sequential_1224896sequential_1224898sequential_1224900sequential_1224902sequential_1224904sequential_1224906sequential_1224908sequential_1224910sequential_1224912sequential_1224914sequential_1224916sequential_1224918sequential_1224920sequential_1224922sequential_1224924sequential_1224926sequential_1224928sequential_1224930sequential_1224932sequential_1224934**
Tin#
!2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
*@
_read_only_resource_inputs"
 	
*0
config_proto 

CPU

GPU2*0J 8? *P
fKRI
G__inference_sequential_layer_call_and_return_conditional_losses_12241862$
"sequential/StatefulPartitionedCall?
$sequential/StatefulPartitionedCall_1StatefulPartitionedCallinputs_1sequential_1224876sequential_1224878sequential_1224880sequential_1224882sequential_1224884sequential_1224886sequential_1224888sequential_1224890sequential_1224892sequential_1224894sequential_1224896sequential_1224898sequential_1224900sequential_1224902sequential_1224904sequential_1224906sequential_1224908sequential_1224910sequential_1224912sequential_1224914sequential_1224916sequential_1224918sequential_1224920sequential_1224922sequential_1224924sequential_1224926sequential_1224928sequential_1224930sequential_1224932sequential_1224934**
Tin#
!2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
*@
_read_only_resource_inputs"
 	
*0
config_proto 

CPU

GPU2*0J 8? *P
fKRI
G__inference_sequential_layer_call_and_return_conditional_losses_12241862&
$sequential/StatefulPartitionedCall_1?
lambda/PartitionedCallPartitionedCall+sequential/StatefulPartitionedCall:output:0-sequential/StatefulPartitionedCall_1:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *L
fGRE
C__inference_lambda_layer_call_and_return_conditional_losses_12244912
lambda/PartitionedCall?
concatenate/PartitionedCallPartitionedCalllambda/PartitionedCall:output:0inputs_2*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *Q
fLRJ
H__inference_concatenate_layer_call_and_return_conditional_losses_12245132
concatenate/PartitionedCall?
dense_2/StatefulPartitionedCallStatefulPartitionedCall$concatenate/PartitionedCall:output:0dense_2_1224970dense_2_1224972*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_2_layer_call_and_return_conditional_losses_12245322!
dense_2/StatefulPartitionedCall?
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_1224928*
_output_shapes
:	?d*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOp?
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?d2!
dense/kernel/Regularizer/Square?
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const?
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/Sum?
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2 
dense/kernel/Regularizer/mul/x?
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/mul?
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_1224932*
_output_shapes

:d
*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp?
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d
2#
!dense_1/kernel/Regularizer/Square?
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const?
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/Sum?
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2"
 dense_1/kernel/Regularizer/mul/x?
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/mul?
IdentityIdentity(dense_2/StatefulPartitionedCall:output:0 ^dense_2/StatefulPartitionedCall#^sequential/StatefulPartitionedCall%^sequential/StatefulPartitionedCall_1*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????:??????????:?????????
::::::::::::::::::::::::::::::::2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2H
"sequential/StatefulPartitionedCall"sequential/StatefulPartitionedCall2L
$sequential/StatefulPartitionedCall_1$sequential/StatefulPartitionedCall_1:\ X
4
_output_shapes"
 :??????????
 
_user_specified_nameinputs:\X
4
_output_shapes"
 :??????????
 
_user_specified_nameinputs:OK
'
_output_shapes
:?????????

 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_1227435

inputs%
!batchnorm_readvariableop_resource)
%batchnorm_mul_readvariableop_resource'
#batchnorm_readvariableop_1_resource'
#batchnorm_readvariableop_2_resource
identity??
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/mul_1?
batchnorm/ReadVariableOp_1ReadVariableOp#batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_1?
batchnorm/mul_2Mul"batchnorm/ReadVariableOp_1:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOp_2ReadVariableOp#batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_2?
batchnorm/subSub"batchnorm/ReadVariableOp_2:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/add_1s
IdentityIdentitybatchnorm/add_1:z:0*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????:::::[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?
?
.__inference_functional_1_layer_call_fn_1223631
input_1
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16

unknown_17

unknown_18

unknown_19

unknown_20

unknown_21

unknown_22

unknown_23

unknown_24
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinput_1unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19
unknown_20
unknown_21
unknown_22
unknown_23
unknown_24*&
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*<
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8? *R
fMRK
I__inference_functional_1_layer_call_and_return_conditional_losses_12235762
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????::::::::::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:] Y
4
_output_shapes"
 :??????????
!
_user_specified_name	input_1
?+
?
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_1227517

inputs
assignmovingavg_1227492
assignmovingavg_1_1227498)
%batchnorm_mul_readvariableop_resource%
!batchnorm_readvariableop_resource
identity??#AssignMovingAvg/AssignSubVariableOp?%AssignMovingAvg_1/AssignSubVariableOp?
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2 
moments/mean/reduction_indices?
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/mean?
moments/StopGradientStopGradientmoments/mean:output:0*
T0**
_output_shapes
:2
moments/StopGradient?
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
moments/SquaredDifference?
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2$
"moments/variance/reduction_indices?
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/variance?
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze?
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze_1?
AssignMovingAvg/decayConst**
_class 
loc:@AssignMovingAvg/1227492*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg/decay?
AssignMovingAvg/ReadVariableOpReadVariableOpassignmovingavg_1227492*
_output_shapes
:*
dtype02 
AssignMovingAvg/ReadVariableOp?
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0**
_class 
loc:@AssignMovingAvg/1227492*
_output_shapes
:2
AssignMovingAvg/sub?
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0**
_class 
loc:@AssignMovingAvg/1227492*
_output_shapes
:2
AssignMovingAvg/mul?
#AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1227492AssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp**
_class 
loc:@AssignMovingAvg/1227492*
_output_shapes
 *
dtype02%
#AssignMovingAvg/AssignSubVariableOp?
AssignMovingAvg_1/decayConst*,
_class"
 loc:@AssignMovingAvg_1/1227498*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg_1/decay?
 AssignMovingAvg_1/ReadVariableOpReadVariableOpassignmovingavg_1_1227498*
_output_shapes
:*
dtype02"
 AssignMovingAvg_1/ReadVariableOp?
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1227498*
_output_shapes
:2
AssignMovingAvg_1/sub?
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1227498*
_output_shapes
:2
AssignMovingAvg_1/mul?
%AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1_1227498AssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*,
_class"
 loc:@AssignMovingAvg_1/1227498*
_output_shapes
 *
dtype02'
%AssignMovingAvg_1/AssignSubVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/mul_1{
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp?
batchnorm/subSub batchnorm/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0$^AssignMovingAvg/AssignSubVariableOp&^AssignMovingAvg_1/AssignSubVariableOp*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????::::2J
#AssignMovingAvg/AssignSubVariableOp#AssignMovingAvg/AssignSubVariableOp2N
%AssignMovingAvg_1/AssignSubVariableOp%AssignMovingAvg_1/AssignSubVariableOp:v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
?
?
.__inference_functional_1_layer_call_fn_1226767

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16

unknown_17

unknown_18

unknown_19

unknown_20

unknown_21

unknown_22

unknown_23

unknown_24
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19
unknown_20
unknown_21
unknown_22
unknown_23
unknown_24*&
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*4
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *R
fMRK
I__inference_functional_1_layer_call_and_return_conditional_losses_12234512
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????::::::::::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:\ X
4
_output_shapes"
 :??????????
 
_user_specified_nameinputs
?
?
7__inference_batch_normalization_1_layer_call_fn_1227277

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_12229672
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????::::22
StatefulPartitionedCallStatefulPartitionedCall:[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?
?
C__inference_conv3d_layer_call_and_return_conditional_losses_1222778

inputs"
conv3d_readvariableop_resource#
biasadd_readvariableop_resource
identity??
Conv3D/ReadVariableOpReadVariableOpconv3d_readvariableop_resource*+
_output_shapes
:?*
dtype02
Conv3D/ReadVariableOp?
Conv3DConv3DinputsConv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
Conv3D?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddConv3D:output:0BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2	
BiasAddp
IdentityIdentityBiasAdd:output:0*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*;
_input_shapes*
(:??????????:::\ X
4
_output_shapes"
 :??????????
 
_user_specified_nameinputs
?

*__inference_conv3d_3_layer_call_fn_1227481

inputs
unknown
	unknown_0
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv3d_3_layer_call_and_return_conditional_losses_12231322
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*:
_input_shapes)
':?????????::22
StatefulPartitionedCallStatefulPartitionedCall:[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
??
?
I__inference_functional_3_layer_call_and_return_conditional_losses_1225520
inputs_0
inputs_1
inputs_2A
=sequential_functional_1_conv3d_conv3d_readvariableop_resourceB
>sequential_functional_1_conv3d_biasadd_readvariableop_resourceG
Csequential_functional_1_batch_normalization_assignmovingavg_1225167I
Esequential_functional_1_batch_normalization_assignmovingavg_1_1225173U
Qsequential_functional_1_batch_normalization_batchnorm_mul_readvariableop_resourceQ
Msequential_functional_1_batch_normalization_batchnorm_readvariableop_resourceC
?sequential_functional_1_conv3d_1_conv3d_readvariableop_resourceD
@sequential_functional_1_conv3d_1_biasadd_readvariableop_resourceI
Esequential_functional_1_batch_normalization_1_assignmovingavg_1225206K
Gsequential_functional_1_batch_normalization_1_assignmovingavg_1_1225212W
Ssequential_functional_1_batch_normalization_1_batchnorm_mul_readvariableop_resourceS
Osequential_functional_1_batch_normalization_1_batchnorm_readvariableop_resourceC
?sequential_functional_1_conv3d_2_conv3d_readvariableop_resourceD
@sequential_functional_1_conv3d_2_biasadd_readvariableop_resourceI
Esequential_functional_1_batch_normalization_2_assignmovingavg_1225245K
Gsequential_functional_1_batch_normalization_2_assignmovingavg_1_1225251W
Ssequential_functional_1_batch_normalization_2_batchnorm_mul_readvariableop_resourceS
Osequential_functional_1_batch_normalization_2_batchnorm_readvariableop_resourceC
?sequential_functional_1_conv3d_3_conv3d_readvariableop_resourceD
@sequential_functional_1_conv3d_3_biasadd_readvariableop_resourceI
Esequential_functional_1_batch_normalization_3_assignmovingavg_1225284K
Gsequential_functional_1_batch_normalization_3_assignmovingavg_1_1225290W
Ssequential_functional_1_batch_normalization_3_batchnorm_mul_readvariableop_resourceS
Osequential_functional_1_batch_normalization_3_batchnorm_readvariableop_resourceA
=sequential_functional_1_layer1_matmul_readvariableop_resourceB
>sequential_functional_1_layer1_biasadd_readvariableop_resource3
/sequential_dense_matmul_readvariableop_resource4
0sequential_dense_biasadd_readvariableop_resource5
1sequential_dense_1_matmul_readvariableop_resource6
2sequential_dense_1_biasadd_readvariableop_resource*
&dense_2_matmul_readvariableop_resource+
'dense_2_biasadd_readvariableop_resource
identity??Osequential/functional_1/batch_normalization/AssignMovingAvg/AssignSubVariableOp?Qsequential/functional_1/batch_normalization/AssignMovingAvg_1/AssignSubVariableOp?Qsequential/functional_1/batch_normalization/AssignMovingAvg_2/AssignSubVariableOp?Qsequential/functional_1/batch_normalization/AssignMovingAvg_3/AssignSubVariableOp?Qsequential/functional_1/batch_normalization_1/AssignMovingAvg/AssignSubVariableOp?Ssequential/functional_1/batch_normalization_1/AssignMovingAvg_1/AssignSubVariableOp?Ssequential/functional_1/batch_normalization_1/AssignMovingAvg_2/AssignSubVariableOp?Ssequential/functional_1/batch_normalization_1/AssignMovingAvg_3/AssignSubVariableOp?Qsequential/functional_1/batch_normalization_2/AssignMovingAvg/AssignSubVariableOp?Ssequential/functional_1/batch_normalization_2/AssignMovingAvg_1/AssignSubVariableOp?Ssequential/functional_1/batch_normalization_2/AssignMovingAvg_2/AssignSubVariableOp?Ssequential/functional_1/batch_normalization_2/AssignMovingAvg_3/AssignSubVariableOp?Qsequential/functional_1/batch_normalization_3/AssignMovingAvg/AssignSubVariableOp?Ssequential/functional_1/batch_normalization_3/AssignMovingAvg_1/AssignSubVariableOp?Ssequential/functional_1/batch_normalization_3/AssignMovingAvg_2/AssignSubVariableOp?Ssequential/functional_1/batch_normalization_3/AssignMovingAvg_3/AssignSubVariableOp?
4sequential/functional_1/conv3d/Conv3D/ReadVariableOpReadVariableOp=sequential_functional_1_conv3d_conv3d_readvariableop_resource*+
_output_shapes
:?*
dtype026
4sequential/functional_1/conv3d/Conv3D/ReadVariableOp?
%sequential/functional_1/conv3d/Conv3DConv3Dinputs_0<sequential/functional_1/conv3d/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2'
%sequential/functional_1/conv3d/Conv3D?
5sequential/functional_1/conv3d/BiasAdd/ReadVariableOpReadVariableOp>sequential_functional_1_conv3d_biasadd_readvariableop_resource*
_output_shapes
:*
dtype027
5sequential/functional_1/conv3d/BiasAdd/ReadVariableOp?
&sequential/functional_1/conv3d/BiasAddBiasAdd.sequential/functional_1/conv3d/Conv3D:output:0=sequential/functional_1/conv3d/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2(
&sequential/functional_1/conv3d/BiasAdd?
Jsequential/functional_1/batch_normalization/moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2L
Jsequential/functional_1/batch_normalization/moments/mean/reduction_indices?
8sequential/functional_1/batch_normalization/moments/meanMean/sequential/functional_1/conv3d/BiasAdd:output:0Ssequential/functional_1/batch_normalization/moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2:
8sequential/functional_1/batch_normalization/moments/mean?
@sequential/functional_1/batch_normalization/moments/StopGradientStopGradientAsequential/functional_1/batch_normalization/moments/mean:output:0*
T0**
_output_shapes
:2B
@sequential/functional_1/batch_normalization/moments/StopGradient?
Esequential/functional_1/batch_normalization/moments/SquaredDifferenceSquaredDifference/sequential/functional_1/conv3d/BiasAdd:output:0Isequential/functional_1/batch_normalization/moments/StopGradient:output:0*
T0*3
_output_shapes!
:?????????2G
Esequential/functional_1/batch_normalization/moments/SquaredDifference?
Nsequential/functional_1/batch_normalization/moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2P
Nsequential/functional_1/batch_normalization/moments/variance/reduction_indices?
<sequential/functional_1/batch_normalization/moments/varianceMeanIsequential/functional_1/batch_normalization/moments/SquaredDifference:z:0Wsequential/functional_1/batch_normalization/moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2>
<sequential/functional_1/batch_normalization/moments/variance?
;sequential/functional_1/batch_normalization/moments/SqueezeSqueezeAsequential/functional_1/batch_normalization/moments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2=
;sequential/functional_1/batch_normalization/moments/Squeeze?
=sequential/functional_1/batch_normalization/moments/Squeeze_1SqueezeEsequential/functional_1/batch_normalization/moments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2?
=sequential/functional_1/batch_normalization/moments/Squeeze_1?
Asequential/functional_1/batch_normalization/AssignMovingAvg/decayConst*V
_classL
JHloc:@sequential/functional_1/batch_normalization/AssignMovingAvg/1225167*
_output_shapes
: *
dtype0*
valueB
 *
?#<2C
Asequential/functional_1/batch_normalization/AssignMovingAvg/decay?
Jsequential/functional_1/batch_normalization/AssignMovingAvg/ReadVariableOpReadVariableOpCsequential_functional_1_batch_normalization_assignmovingavg_1225167*
_output_shapes
:*
dtype02L
Jsequential/functional_1/batch_normalization/AssignMovingAvg/ReadVariableOp?
?sequential/functional_1/batch_normalization/AssignMovingAvg/subSubRsequential/functional_1/batch_normalization/AssignMovingAvg/ReadVariableOp:value:0Dsequential/functional_1/batch_normalization/moments/Squeeze:output:0*
T0*V
_classL
JHloc:@sequential/functional_1/batch_normalization/AssignMovingAvg/1225167*
_output_shapes
:2A
?sequential/functional_1/batch_normalization/AssignMovingAvg/sub?
?sequential/functional_1/batch_normalization/AssignMovingAvg/mulMulCsequential/functional_1/batch_normalization/AssignMovingAvg/sub:z:0Jsequential/functional_1/batch_normalization/AssignMovingAvg/decay:output:0*
T0*V
_classL
JHloc:@sequential/functional_1/batch_normalization/AssignMovingAvg/1225167*
_output_shapes
:2A
?sequential/functional_1/batch_normalization/AssignMovingAvg/mul?
Osequential/functional_1/batch_normalization/AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpCsequential_functional_1_batch_normalization_assignmovingavg_1225167Csequential/functional_1/batch_normalization/AssignMovingAvg/mul:z:0K^sequential/functional_1/batch_normalization/AssignMovingAvg/ReadVariableOp*V
_classL
JHloc:@sequential/functional_1/batch_normalization/AssignMovingAvg/1225167*
_output_shapes
 *
dtype02Q
Osequential/functional_1/batch_normalization/AssignMovingAvg/AssignSubVariableOp?
Csequential/functional_1/batch_normalization/AssignMovingAvg_1/decayConst*X
_classN
LJloc:@sequential/functional_1/batch_normalization/AssignMovingAvg_1/1225173*
_output_shapes
: *
dtype0*
valueB
 *
?#<2E
Csequential/functional_1/batch_normalization/AssignMovingAvg_1/decay?
Lsequential/functional_1/batch_normalization/AssignMovingAvg_1/ReadVariableOpReadVariableOpEsequential_functional_1_batch_normalization_assignmovingavg_1_1225173*
_output_shapes
:*
dtype02N
Lsequential/functional_1/batch_normalization/AssignMovingAvg_1/ReadVariableOp?
Asequential/functional_1/batch_normalization/AssignMovingAvg_1/subSubTsequential/functional_1/batch_normalization/AssignMovingAvg_1/ReadVariableOp:value:0Fsequential/functional_1/batch_normalization/moments/Squeeze_1:output:0*
T0*X
_classN
LJloc:@sequential/functional_1/batch_normalization/AssignMovingAvg_1/1225173*
_output_shapes
:2C
Asequential/functional_1/batch_normalization/AssignMovingAvg_1/sub?
Asequential/functional_1/batch_normalization/AssignMovingAvg_1/mulMulEsequential/functional_1/batch_normalization/AssignMovingAvg_1/sub:z:0Lsequential/functional_1/batch_normalization/AssignMovingAvg_1/decay:output:0*
T0*X
_classN
LJloc:@sequential/functional_1/batch_normalization/AssignMovingAvg_1/1225173*
_output_shapes
:2C
Asequential/functional_1/batch_normalization/AssignMovingAvg_1/mul?
Qsequential/functional_1/batch_normalization/AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOpEsequential_functional_1_batch_normalization_assignmovingavg_1_1225173Esequential/functional_1/batch_normalization/AssignMovingAvg_1/mul:z:0M^sequential/functional_1/batch_normalization/AssignMovingAvg_1/ReadVariableOp*X
_classN
LJloc:@sequential/functional_1/batch_normalization/AssignMovingAvg_1/1225173*
_output_shapes
 *
dtype02S
Qsequential/functional_1/batch_normalization/AssignMovingAvg_1/AssignSubVariableOp?
;sequential/functional_1/batch_normalization/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2=
;sequential/functional_1/batch_normalization/batchnorm/add/y?
9sequential/functional_1/batch_normalization/batchnorm/addAddV2Fsequential/functional_1/batch_normalization/moments/Squeeze_1:output:0Dsequential/functional_1/batch_normalization/batchnorm/add/y:output:0*
T0*
_output_shapes
:2;
9sequential/functional_1/batch_normalization/batchnorm/add?
;sequential/functional_1/batch_normalization/batchnorm/RsqrtRsqrt=sequential/functional_1/batch_normalization/batchnorm/add:z:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization/batchnorm/Rsqrt?
Hsequential/functional_1/batch_normalization/batchnorm/mul/ReadVariableOpReadVariableOpQsequential_functional_1_batch_normalization_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02J
Hsequential/functional_1/batch_normalization/batchnorm/mul/ReadVariableOp?
9sequential/functional_1/batch_normalization/batchnorm/mulMul?sequential/functional_1/batch_normalization/batchnorm/Rsqrt:y:0Psequential/functional_1/batch_normalization/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2;
9sequential/functional_1/batch_normalization/batchnorm/mul?
;sequential/functional_1/batch_normalization/batchnorm/mul_1Mul/sequential/functional_1/conv3d/BiasAdd:output:0=sequential/functional_1/batch_normalization/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2=
;sequential/functional_1/batch_normalization/batchnorm/mul_1?
;sequential/functional_1/batch_normalization/batchnorm/mul_2MulDsequential/functional_1/batch_normalization/moments/Squeeze:output:0=sequential/functional_1/batch_normalization/batchnorm/mul:z:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization/batchnorm/mul_2?
Dsequential/functional_1/batch_normalization/batchnorm/ReadVariableOpReadVariableOpMsequential_functional_1_batch_normalization_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02F
Dsequential/functional_1/batch_normalization/batchnorm/ReadVariableOp?
9sequential/functional_1/batch_normalization/batchnorm/subSubLsequential/functional_1/batch_normalization/batchnorm/ReadVariableOp:value:0?sequential/functional_1/batch_normalization/batchnorm/mul_2:z:0*
T0*
_output_shapes
:2;
9sequential/functional_1/batch_normalization/batchnorm/sub?
;sequential/functional_1/batch_normalization/batchnorm/add_1AddV2?sequential/functional_1/batch_normalization/batchnorm/mul_1:z:0=sequential/functional_1/batch_normalization/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2=
;sequential/functional_1/batch_normalization/batchnorm/add_1?
6sequential/functional_1/conv3d_1/Conv3D/ReadVariableOpReadVariableOp?sequential_functional_1_conv3d_1_conv3d_readvariableop_resource**
_output_shapes
:*
dtype028
6sequential/functional_1/conv3d_1/Conv3D/ReadVariableOp?
'sequential/functional_1/conv3d_1/Conv3DConv3D?sequential/functional_1/batch_normalization/batchnorm/add_1:z:0>sequential/functional_1/conv3d_1/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2)
'sequential/functional_1/conv3d_1/Conv3D?
7sequential/functional_1/conv3d_1/BiasAdd/ReadVariableOpReadVariableOp@sequential_functional_1_conv3d_1_biasadd_readvariableop_resource*
_output_shapes
:*
dtype029
7sequential/functional_1/conv3d_1/BiasAdd/ReadVariableOp?
(sequential/functional_1/conv3d_1/BiasAddBiasAdd0sequential/functional_1/conv3d_1/Conv3D:output:0?sequential/functional_1/conv3d_1/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2*
(sequential/functional_1/conv3d_1/BiasAdd?
$sequential/functional_1/conv3d_1/EluElu1sequential/functional_1/conv3d_1/BiasAdd:output:0*
T0*3
_output_shapes!
:?????????2&
$sequential/functional_1/conv3d_1/Elu?
Lsequential/functional_1/batch_normalization_1/moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2N
Lsequential/functional_1/batch_normalization_1/moments/mean/reduction_indices?
:sequential/functional_1/batch_normalization_1/moments/meanMean2sequential/functional_1/conv3d_1/Elu:activations:0Usequential/functional_1/batch_normalization_1/moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2<
:sequential/functional_1/batch_normalization_1/moments/mean?
Bsequential/functional_1/batch_normalization_1/moments/StopGradientStopGradientCsequential/functional_1/batch_normalization_1/moments/mean:output:0*
T0**
_output_shapes
:2D
Bsequential/functional_1/batch_normalization_1/moments/StopGradient?
Gsequential/functional_1/batch_normalization_1/moments/SquaredDifferenceSquaredDifference2sequential/functional_1/conv3d_1/Elu:activations:0Ksequential/functional_1/batch_normalization_1/moments/StopGradient:output:0*
T0*3
_output_shapes!
:?????????2I
Gsequential/functional_1/batch_normalization_1/moments/SquaredDifference?
Psequential/functional_1/batch_normalization_1/moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2R
Psequential/functional_1/batch_normalization_1/moments/variance/reduction_indices?
>sequential/functional_1/batch_normalization_1/moments/varianceMeanKsequential/functional_1/batch_normalization_1/moments/SquaredDifference:z:0Ysequential/functional_1/batch_normalization_1/moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2@
>sequential/functional_1/batch_normalization_1/moments/variance?
=sequential/functional_1/batch_normalization_1/moments/SqueezeSqueezeCsequential/functional_1/batch_normalization_1/moments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2?
=sequential/functional_1/batch_normalization_1/moments/Squeeze?
?sequential/functional_1/batch_normalization_1/moments/Squeeze_1SqueezeGsequential/functional_1/batch_normalization_1/moments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2A
?sequential/functional_1/batch_normalization_1/moments/Squeeze_1?
Csequential/functional_1/batch_normalization_1/AssignMovingAvg/decayConst*X
_classN
LJloc:@sequential/functional_1/batch_normalization_1/AssignMovingAvg/1225206*
_output_shapes
: *
dtype0*
valueB
 *
?#<2E
Csequential/functional_1/batch_normalization_1/AssignMovingAvg/decay?
Lsequential/functional_1/batch_normalization_1/AssignMovingAvg/ReadVariableOpReadVariableOpEsequential_functional_1_batch_normalization_1_assignmovingavg_1225206*
_output_shapes
:*
dtype02N
Lsequential/functional_1/batch_normalization_1/AssignMovingAvg/ReadVariableOp?
Asequential/functional_1/batch_normalization_1/AssignMovingAvg/subSubTsequential/functional_1/batch_normalization_1/AssignMovingAvg/ReadVariableOp:value:0Fsequential/functional_1/batch_normalization_1/moments/Squeeze:output:0*
T0*X
_classN
LJloc:@sequential/functional_1/batch_normalization_1/AssignMovingAvg/1225206*
_output_shapes
:2C
Asequential/functional_1/batch_normalization_1/AssignMovingAvg/sub?
Asequential/functional_1/batch_normalization_1/AssignMovingAvg/mulMulEsequential/functional_1/batch_normalization_1/AssignMovingAvg/sub:z:0Lsequential/functional_1/batch_normalization_1/AssignMovingAvg/decay:output:0*
T0*X
_classN
LJloc:@sequential/functional_1/batch_normalization_1/AssignMovingAvg/1225206*
_output_shapes
:2C
Asequential/functional_1/batch_normalization_1/AssignMovingAvg/mul?
Qsequential/functional_1/batch_normalization_1/AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpEsequential_functional_1_batch_normalization_1_assignmovingavg_1225206Esequential/functional_1/batch_normalization_1/AssignMovingAvg/mul:z:0M^sequential/functional_1/batch_normalization_1/AssignMovingAvg/ReadVariableOp*X
_classN
LJloc:@sequential/functional_1/batch_normalization_1/AssignMovingAvg/1225206*
_output_shapes
 *
dtype02S
Qsequential/functional_1/batch_normalization_1/AssignMovingAvg/AssignSubVariableOp?
Esequential/functional_1/batch_normalization_1/AssignMovingAvg_1/decayConst*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_1/AssignMovingAvg_1/1225212*
_output_shapes
: *
dtype0*
valueB
 *
?#<2G
Esequential/functional_1/batch_normalization_1/AssignMovingAvg_1/decay?
Nsequential/functional_1/batch_normalization_1/AssignMovingAvg_1/ReadVariableOpReadVariableOpGsequential_functional_1_batch_normalization_1_assignmovingavg_1_1225212*
_output_shapes
:*
dtype02P
Nsequential/functional_1/batch_normalization_1/AssignMovingAvg_1/ReadVariableOp?
Csequential/functional_1/batch_normalization_1/AssignMovingAvg_1/subSubVsequential/functional_1/batch_normalization_1/AssignMovingAvg_1/ReadVariableOp:value:0Hsequential/functional_1/batch_normalization_1/moments/Squeeze_1:output:0*
T0*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_1/AssignMovingAvg_1/1225212*
_output_shapes
:2E
Csequential/functional_1/batch_normalization_1/AssignMovingAvg_1/sub?
Csequential/functional_1/batch_normalization_1/AssignMovingAvg_1/mulMulGsequential/functional_1/batch_normalization_1/AssignMovingAvg_1/sub:z:0Nsequential/functional_1/batch_normalization_1/AssignMovingAvg_1/decay:output:0*
T0*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_1/AssignMovingAvg_1/1225212*
_output_shapes
:2E
Csequential/functional_1/batch_normalization_1/AssignMovingAvg_1/mul?
Ssequential/functional_1/batch_normalization_1/AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOpGsequential_functional_1_batch_normalization_1_assignmovingavg_1_1225212Gsequential/functional_1/batch_normalization_1/AssignMovingAvg_1/mul:z:0O^sequential/functional_1/batch_normalization_1/AssignMovingAvg_1/ReadVariableOp*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_1/AssignMovingAvg_1/1225212*
_output_shapes
 *
dtype02U
Ssequential/functional_1/batch_normalization_1/AssignMovingAvg_1/AssignSubVariableOp?
=sequential/functional_1/batch_normalization_1/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2?
=sequential/functional_1/batch_normalization_1/batchnorm/add/y?
;sequential/functional_1/batch_normalization_1/batchnorm/addAddV2Hsequential/functional_1/batch_normalization_1/moments/Squeeze_1:output:0Fsequential/functional_1/batch_normalization_1/batchnorm/add/y:output:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization_1/batchnorm/add?
=sequential/functional_1/batch_normalization_1/batchnorm/RsqrtRsqrt?sequential/functional_1/batch_normalization_1/batchnorm/add:z:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_1/batchnorm/Rsqrt?
Jsequential/functional_1/batch_normalization_1/batchnorm/mul/ReadVariableOpReadVariableOpSsequential_functional_1_batch_normalization_1_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02L
Jsequential/functional_1/batch_normalization_1/batchnorm/mul/ReadVariableOp?
;sequential/functional_1/batch_normalization_1/batchnorm/mulMulAsequential/functional_1/batch_normalization_1/batchnorm/Rsqrt:y:0Rsequential/functional_1/batch_normalization_1/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization_1/batchnorm/mul?
=sequential/functional_1/batch_normalization_1/batchnorm/mul_1Mul2sequential/functional_1/conv3d_1/Elu:activations:0?sequential/functional_1/batch_normalization_1/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2?
=sequential/functional_1/batch_normalization_1/batchnorm/mul_1?
=sequential/functional_1/batch_normalization_1/batchnorm/mul_2MulFsequential/functional_1/batch_normalization_1/moments/Squeeze:output:0?sequential/functional_1/batch_normalization_1/batchnorm/mul:z:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_1/batchnorm/mul_2?
Fsequential/functional_1/batch_normalization_1/batchnorm/ReadVariableOpReadVariableOpOsequential_functional_1_batch_normalization_1_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02H
Fsequential/functional_1/batch_normalization_1/batchnorm/ReadVariableOp?
;sequential/functional_1/batch_normalization_1/batchnorm/subSubNsequential/functional_1/batch_normalization_1/batchnorm/ReadVariableOp:value:0Asequential/functional_1/batch_normalization_1/batchnorm/mul_2:z:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization_1/batchnorm/sub?
=sequential/functional_1/batch_normalization_1/batchnorm/add_1AddV2Asequential/functional_1/batch_normalization_1/batchnorm/mul_1:z:0?sequential/functional_1/batch_normalization_1/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2?
=sequential/functional_1/batch_normalization_1/batchnorm/add_1?
6sequential/functional_1/conv3d_2/Conv3D/ReadVariableOpReadVariableOp?sequential_functional_1_conv3d_2_conv3d_readvariableop_resource**
_output_shapes
:*
dtype028
6sequential/functional_1/conv3d_2/Conv3D/ReadVariableOp?
'sequential/functional_1/conv3d_2/Conv3DConv3DAsequential/functional_1/batch_normalization_1/batchnorm/add_1:z:0>sequential/functional_1/conv3d_2/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2)
'sequential/functional_1/conv3d_2/Conv3D?
7sequential/functional_1/conv3d_2/BiasAdd/ReadVariableOpReadVariableOp@sequential_functional_1_conv3d_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype029
7sequential/functional_1/conv3d_2/BiasAdd/ReadVariableOp?
(sequential/functional_1/conv3d_2/BiasAddBiasAdd0sequential/functional_1/conv3d_2/Conv3D:output:0?sequential/functional_1/conv3d_2/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2*
(sequential/functional_1/conv3d_2/BiasAdd?
$sequential/functional_1/conv3d_2/EluElu1sequential/functional_1/conv3d_2/BiasAdd:output:0*
T0*3
_output_shapes!
:?????????2&
$sequential/functional_1/conv3d_2/Elu?
Lsequential/functional_1/batch_normalization_2/moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2N
Lsequential/functional_1/batch_normalization_2/moments/mean/reduction_indices?
:sequential/functional_1/batch_normalization_2/moments/meanMean2sequential/functional_1/conv3d_2/Elu:activations:0Usequential/functional_1/batch_normalization_2/moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2<
:sequential/functional_1/batch_normalization_2/moments/mean?
Bsequential/functional_1/batch_normalization_2/moments/StopGradientStopGradientCsequential/functional_1/batch_normalization_2/moments/mean:output:0*
T0**
_output_shapes
:2D
Bsequential/functional_1/batch_normalization_2/moments/StopGradient?
Gsequential/functional_1/batch_normalization_2/moments/SquaredDifferenceSquaredDifference2sequential/functional_1/conv3d_2/Elu:activations:0Ksequential/functional_1/batch_normalization_2/moments/StopGradient:output:0*
T0*3
_output_shapes!
:?????????2I
Gsequential/functional_1/batch_normalization_2/moments/SquaredDifference?
Psequential/functional_1/batch_normalization_2/moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2R
Psequential/functional_1/batch_normalization_2/moments/variance/reduction_indices?
>sequential/functional_1/batch_normalization_2/moments/varianceMeanKsequential/functional_1/batch_normalization_2/moments/SquaredDifference:z:0Ysequential/functional_1/batch_normalization_2/moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2@
>sequential/functional_1/batch_normalization_2/moments/variance?
=sequential/functional_1/batch_normalization_2/moments/SqueezeSqueezeCsequential/functional_1/batch_normalization_2/moments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2?
=sequential/functional_1/batch_normalization_2/moments/Squeeze?
?sequential/functional_1/batch_normalization_2/moments/Squeeze_1SqueezeGsequential/functional_1/batch_normalization_2/moments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2A
?sequential/functional_1/batch_normalization_2/moments/Squeeze_1?
Csequential/functional_1/batch_normalization_2/AssignMovingAvg/decayConst*X
_classN
LJloc:@sequential/functional_1/batch_normalization_2/AssignMovingAvg/1225245*
_output_shapes
: *
dtype0*
valueB
 *
?#<2E
Csequential/functional_1/batch_normalization_2/AssignMovingAvg/decay?
Lsequential/functional_1/batch_normalization_2/AssignMovingAvg/ReadVariableOpReadVariableOpEsequential_functional_1_batch_normalization_2_assignmovingavg_1225245*
_output_shapes
:*
dtype02N
Lsequential/functional_1/batch_normalization_2/AssignMovingAvg/ReadVariableOp?
Asequential/functional_1/batch_normalization_2/AssignMovingAvg/subSubTsequential/functional_1/batch_normalization_2/AssignMovingAvg/ReadVariableOp:value:0Fsequential/functional_1/batch_normalization_2/moments/Squeeze:output:0*
T0*X
_classN
LJloc:@sequential/functional_1/batch_normalization_2/AssignMovingAvg/1225245*
_output_shapes
:2C
Asequential/functional_1/batch_normalization_2/AssignMovingAvg/sub?
Asequential/functional_1/batch_normalization_2/AssignMovingAvg/mulMulEsequential/functional_1/batch_normalization_2/AssignMovingAvg/sub:z:0Lsequential/functional_1/batch_normalization_2/AssignMovingAvg/decay:output:0*
T0*X
_classN
LJloc:@sequential/functional_1/batch_normalization_2/AssignMovingAvg/1225245*
_output_shapes
:2C
Asequential/functional_1/batch_normalization_2/AssignMovingAvg/mul?
Qsequential/functional_1/batch_normalization_2/AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpEsequential_functional_1_batch_normalization_2_assignmovingavg_1225245Esequential/functional_1/batch_normalization_2/AssignMovingAvg/mul:z:0M^sequential/functional_1/batch_normalization_2/AssignMovingAvg/ReadVariableOp*X
_classN
LJloc:@sequential/functional_1/batch_normalization_2/AssignMovingAvg/1225245*
_output_shapes
 *
dtype02S
Qsequential/functional_1/batch_normalization_2/AssignMovingAvg/AssignSubVariableOp?
Esequential/functional_1/batch_normalization_2/AssignMovingAvg_1/decayConst*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_2/AssignMovingAvg_1/1225251*
_output_shapes
: *
dtype0*
valueB
 *
?#<2G
Esequential/functional_1/batch_normalization_2/AssignMovingAvg_1/decay?
Nsequential/functional_1/batch_normalization_2/AssignMovingAvg_1/ReadVariableOpReadVariableOpGsequential_functional_1_batch_normalization_2_assignmovingavg_1_1225251*
_output_shapes
:*
dtype02P
Nsequential/functional_1/batch_normalization_2/AssignMovingAvg_1/ReadVariableOp?
Csequential/functional_1/batch_normalization_2/AssignMovingAvg_1/subSubVsequential/functional_1/batch_normalization_2/AssignMovingAvg_1/ReadVariableOp:value:0Hsequential/functional_1/batch_normalization_2/moments/Squeeze_1:output:0*
T0*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_2/AssignMovingAvg_1/1225251*
_output_shapes
:2E
Csequential/functional_1/batch_normalization_2/AssignMovingAvg_1/sub?
Csequential/functional_1/batch_normalization_2/AssignMovingAvg_1/mulMulGsequential/functional_1/batch_normalization_2/AssignMovingAvg_1/sub:z:0Nsequential/functional_1/batch_normalization_2/AssignMovingAvg_1/decay:output:0*
T0*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_2/AssignMovingAvg_1/1225251*
_output_shapes
:2E
Csequential/functional_1/batch_normalization_2/AssignMovingAvg_1/mul?
Ssequential/functional_1/batch_normalization_2/AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOpGsequential_functional_1_batch_normalization_2_assignmovingavg_1_1225251Gsequential/functional_1/batch_normalization_2/AssignMovingAvg_1/mul:z:0O^sequential/functional_1/batch_normalization_2/AssignMovingAvg_1/ReadVariableOp*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_2/AssignMovingAvg_1/1225251*
_output_shapes
 *
dtype02U
Ssequential/functional_1/batch_normalization_2/AssignMovingAvg_1/AssignSubVariableOp?
=sequential/functional_1/batch_normalization_2/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2?
=sequential/functional_1/batch_normalization_2/batchnorm/add/y?
;sequential/functional_1/batch_normalization_2/batchnorm/addAddV2Hsequential/functional_1/batch_normalization_2/moments/Squeeze_1:output:0Fsequential/functional_1/batch_normalization_2/batchnorm/add/y:output:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization_2/batchnorm/add?
=sequential/functional_1/batch_normalization_2/batchnorm/RsqrtRsqrt?sequential/functional_1/batch_normalization_2/batchnorm/add:z:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_2/batchnorm/Rsqrt?
Jsequential/functional_1/batch_normalization_2/batchnorm/mul/ReadVariableOpReadVariableOpSsequential_functional_1_batch_normalization_2_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02L
Jsequential/functional_1/batch_normalization_2/batchnorm/mul/ReadVariableOp?
;sequential/functional_1/batch_normalization_2/batchnorm/mulMulAsequential/functional_1/batch_normalization_2/batchnorm/Rsqrt:y:0Rsequential/functional_1/batch_normalization_2/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization_2/batchnorm/mul?
=sequential/functional_1/batch_normalization_2/batchnorm/mul_1Mul2sequential/functional_1/conv3d_2/Elu:activations:0?sequential/functional_1/batch_normalization_2/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2?
=sequential/functional_1/batch_normalization_2/batchnorm/mul_1?
=sequential/functional_1/batch_normalization_2/batchnorm/mul_2MulFsequential/functional_1/batch_normalization_2/moments/Squeeze:output:0?sequential/functional_1/batch_normalization_2/batchnorm/mul:z:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_2/batchnorm/mul_2?
Fsequential/functional_1/batch_normalization_2/batchnorm/ReadVariableOpReadVariableOpOsequential_functional_1_batch_normalization_2_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02H
Fsequential/functional_1/batch_normalization_2/batchnorm/ReadVariableOp?
;sequential/functional_1/batch_normalization_2/batchnorm/subSubNsequential/functional_1/batch_normalization_2/batchnorm/ReadVariableOp:value:0Asequential/functional_1/batch_normalization_2/batchnorm/mul_2:z:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization_2/batchnorm/sub?
=sequential/functional_1/batch_normalization_2/batchnorm/add_1AddV2Asequential/functional_1/batch_normalization_2/batchnorm/mul_1:z:0?sequential/functional_1/batch_normalization_2/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2?
=sequential/functional_1/batch_normalization_2/batchnorm/add_1?
6sequential/functional_1/conv3d_3/Conv3D/ReadVariableOpReadVariableOp?sequential_functional_1_conv3d_3_conv3d_readvariableop_resource**
_output_shapes
:*
dtype028
6sequential/functional_1/conv3d_3/Conv3D/ReadVariableOp?
'sequential/functional_1/conv3d_3/Conv3DConv3DAsequential/functional_1/batch_normalization_2/batchnorm/add_1:z:0>sequential/functional_1/conv3d_3/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2)
'sequential/functional_1/conv3d_3/Conv3D?
7sequential/functional_1/conv3d_3/BiasAdd/ReadVariableOpReadVariableOp@sequential_functional_1_conv3d_3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype029
7sequential/functional_1/conv3d_3/BiasAdd/ReadVariableOp?
(sequential/functional_1/conv3d_3/BiasAddBiasAdd0sequential/functional_1/conv3d_3/Conv3D:output:0?sequential/functional_1/conv3d_3/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2*
(sequential/functional_1/conv3d_3/BiasAdd?
$sequential/functional_1/conv3d_3/EluElu1sequential/functional_1/conv3d_3/BiasAdd:output:0*
T0*3
_output_shapes!
:?????????2&
$sequential/functional_1/conv3d_3/Elu?
Lsequential/functional_1/batch_normalization_3/moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2N
Lsequential/functional_1/batch_normalization_3/moments/mean/reduction_indices?
:sequential/functional_1/batch_normalization_3/moments/meanMean2sequential/functional_1/conv3d_3/Elu:activations:0Usequential/functional_1/batch_normalization_3/moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2<
:sequential/functional_1/batch_normalization_3/moments/mean?
Bsequential/functional_1/batch_normalization_3/moments/StopGradientStopGradientCsequential/functional_1/batch_normalization_3/moments/mean:output:0*
T0**
_output_shapes
:2D
Bsequential/functional_1/batch_normalization_3/moments/StopGradient?
Gsequential/functional_1/batch_normalization_3/moments/SquaredDifferenceSquaredDifference2sequential/functional_1/conv3d_3/Elu:activations:0Ksequential/functional_1/batch_normalization_3/moments/StopGradient:output:0*
T0*3
_output_shapes!
:?????????2I
Gsequential/functional_1/batch_normalization_3/moments/SquaredDifference?
Psequential/functional_1/batch_normalization_3/moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2R
Psequential/functional_1/batch_normalization_3/moments/variance/reduction_indices?
>sequential/functional_1/batch_normalization_3/moments/varianceMeanKsequential/functional_1/batch_normalization_3/moments/SquaredDifference:z:0Ysequential/functional_1/batch_normalization_3/moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2@
>sequential/functional_1/batch_normalization_3/moments/variance?
=sequential/functional_1/batch_normalization_3/moments/SqueezeSqueezeCsequential/functional_1/batch_normalization_3/moments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2?
=sequential/functional_1/batch_normalization_3/moments/Squeeze?
?sequential/functional_1/batch_normalization_3/moments/Squeeze_1SqueezeGsequential/functional_1/batch_normalization_3/moments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2A
?sequential/functional_1/batch_normalization_3/moments/Squeeze_1?
Csequential/functional_1/batch_normalization_3/AssignMovingAvg/decayConst*X
_classN
LJloc:@sequential/functional_1/batch_normalization_3/AssignMovingAvg/1225284*
_output_shapes
: *
dtype0*
valueB
 *
?#<2E
Csequential/functional_1/batch_normalization_3/AssignMovingAvg/decay?
Lsequential/functional_1/batch_normalization_3/AssignMovingAvg/ReadVariableOpReadVariableOpEsequential_functional_1_batch_normalization_3_assignmovingavg_1225284*
_output_shapes
:*
dtype02N
Lsequential/functional_1/batch_normalization_3/AssignMovingAvg/ReadVariableOp?
Asequential/functional_1/batch_normalization_3/AssignMovingAvg/subSubTsequential/functional_1/batch_normalization_3/AssignMovingAvg/ReadVariableOp:value:0Fsequential/functional_1/batch_normalization_3/moments/Squeeze:output:0*
T0*X
_classN
LJloc:@sequential/functional_1/batch_normalization_3/AssignMovingAvg/1225284*
_output_shapes
:2C
Asequential/functional_1/batch_normalization_3/AssignMovingAvg/sub?
Asequential/functional_1/batch_normalization_3/AssignMovingAvg/mulMulEsequential/functional_1/batch_normalization_3/AssignMovingAvg/sub:z:0Lsequential/functional_1/batch_normalization_3/AssignMovingAvg/decay:output:0*
T0*X
_classN
LJloc:@sequential/functional_1/batch_normalization_3/AssignMovingAvg/1225284*
_output_shapes
:2C
Asequential/functional_1/batch_normalization_3/AssignMovingAvg/mul?
Qsequential/functional_1/batch_normalization_3/AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpEsequential_functional_1_batch_normalization_3_assignmovingavg_1225284Esequential/functional_1/batch_normalization_3/AssignMovingAvg/mul:z:0M^sequential/functional_1/batch_normalization_3/AssignMovingAvg/ReadVariableOp*X
_classN
LJloc:@sequential/functional_1/batch_normalization_3/AssignMovingAvg/1225284*
_output_shapes
 *
dtype02S
Qsequential/functional_1/batch_normalization_3/AssignMovingAvg/AssignSubVariableOp?
Esequential/functional_1/batch_normalization_3/AssignMovingAvg_1/decayConst*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_3/AssignMovingAvg_1/1225290*
_output_shapes
: *
dtype0*
valueB
 *
?#<2G
Esequential/functional_1/batch_normalization_3/AssignMovingAvg_1/decay?
Nsequential/functional_1/batch_normalization_3/AssignMovingAvg_1/ReadVariableOpReadVariableOpGsequential_functional_1_batch_normalization_3_assignmovingavg_1_1225290*
_output_shapes
:*
dtype02P
Nsequential/functional_1/batch_normalization_3/AssignMovingAvg_1/ReadVariableOp?
Csequential/functional_1/batch_normalization_3/AssignMovingAvg_1/subSubVsequential/functional_1/batch_normalization_3/AssignMovingAvg_1/ReadVariableOp:value:0Hsequential/functional_1/batch_normalization_3/moments/Squeeze_1:output:0*
T0*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_3/AssignMovingAvg_1/1225290*
_output_shapes
:2E
Csequential/functional_1/batch_normalization_3/AssignMovingAvg_1/sub?
Csequential/functional_1/batch_normalization_3/AssignMovingAvg_1/mulMulGsequential/functional_1/batch_normalization_3/AssignMovingAvg_1/sub:z:0Nsequential/functional_1/batch_normalization_3/AssignMovingAvg_1/decay:output:0*
T0*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_3/AssignMovingAvg_1/1225290*
_output_shapes
:2E
Csequential/functional_1/batch_normalization_3/AssignMovingAvg_1/mul?
Ssequential/functional_1/batch_normalization_3/AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOpGsequential_functional_1_batch_normalization_3_assignmovingavg_1_1225290Gsequential/functional_1/batch_normalization_3/AssignMovingAvg_1/mul:z:0O^sequential/functional_1/batch_normalization_3/AssignMovingAvg_1/ReadVariableOp*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_3/AssignMovingAvg_1/1225290*
_output_shapes
 *
dtype02U
Ssequential/functional_1/batch_normalization_3/AssignMovingAvg_1/AssignSubVariableOp?
=sequential/functional_1/batch_normalization_3/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2?
=sequential/functional_1/batch_normalization_3/batchnorm/add/y?
;sequential/functional_1/batch_normalization_3/batchnorm/addAddV2Hsequential/functional_1/batch_normalization_3/moments/Squeeze_1:output:0Fsequential/functional_1/batch_normalization_3/batchnorm/add/y:output:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization_3/batchnorm/add?
=sequential/functional_1/batch_normalization_3/batchnorm/RsqrtRsqrt?sequential/functional_1/batch_normalization_3/batchnorm/add:z:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_3/batchnorm/Rsqrt?
Jsequential/functional_1/batch_normalization_3/batchnorm/mul/ReadVariableOpReadVariableOpSsequential_functional_1_batch_normalization_3_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02L
Jsequential/functional_1/batch_normalization_3/batchnorm/mul/ReadVariableOp?
;sequential/functional_1/batch_normalization_3/batchnorm/mulMulAsequential/functional_1/batch_normalization_3/batchnorm/Rsqrt:y:0Rsequential/functional_1/batch_normalization_3/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization_3/batchnorm/mul?
=sequential/functional_1/batch_normalization_3/batchnorm/mul_1Mul2sequential/functional_1/conv3d_3/Elu:activations:0?sequential/functional_1/batch_normalization_3/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2?
=sequential/functional_1/batch_normalization_3/batchnorm/mul_1?
=sequential/functional_1/batch_normalization_3/batchnorm/mul_2MulFsequential/functional_1/batch_normalization_3/moments/Squeeze:output:0?sequential/functional_1/batch_normalization_3/batchnorm/mul:z:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_3/batchnorm/mul_2?
Fsequential/functional_1/batch_normalization_3/batchnorm/ReadVariableOpReadVariableOpOsequential_functional_1_batch_normalization_3_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02H
Fsequential/functional_1/batch_normalization_3/batchnorm/ReadVariableOp?
;sequential/functional_1/batch_normalization_3/batchnorm/subSubNsequential/functional_1/batch_normalization_3/batchnorm/ReadVariableOp:value:0Asequential/functional_1/batch_normalization_3/batchnorm/mul_2:z:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization_3/batchnorm/sub?
=sequential/functional_1/batch_normalization_3/batchnorm/add_1AddV2Asequential/functional_1/batch_normalization_3/batchnorm/mul_1:z:0?sequential/functional_1/batch_normalization_3/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2?
=sequential/functional_1/batch_normalization_3/batchnorm/add_1?
3sequential/functional_1/average_pooling3d/AvgPool3D	AvgPool3DAsequential/functional_1/batch_normalization_3/batchnorm/add_1:z:0*
T0*3
_output_shapes!
:?????????*
ksize	
*
paddingVALID*
strides	
25
3sequential/functional_1/average_pooling3d/AvgPool3D?
%sequential/functional_1/flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"????   2'
%sequential/functional_1/flatten/Const?
'sequential/functional_1/flatten/ReshapeReshape<sequential/functional_1/average_pooling3d/AvgPool3D:output:0.sequential/functional_1/flatten/Const:output:0*
T0*(
_output_shapes
:??????????
2)
'sequential/functional_1/flatten/Reshape?
-sequential/functional_1/dropout/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *UU??2/
-sequential/functional_1/dropout/dropout/Const?
+sequential/functional_1/dropout/dropout/MulMul0sequential/functional_1/flatten/Reshape:output:06sequential/functional_1/dropout/dropout/Const:output:0*
T0*(
_output_shapes
:??????????
2-
+sequential/functional_1/dropout/dropout/Mul?
-sequential/functional_1/dropout/dropout/ShapeShape0sequential/functional_1/flatten/Reshape:output:0*
T0*
_output_shapes
:2/
-sequential/functional_1/dropout/dropout/Shape?
Dsequential/functional_1/dropout/dropout/random_uniform/RandomUniformRandomUniform6sequential/functional_1/dropout/dropout/Shape:output:0*
T0*(
_output_shapes
:??????????
*
dtype02F
Dsequential/functional_1/dropout/dropout/random_uniform/RandomUniform?
6sequential/functional_1/dropout/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>28
6sequential/functional_1/dropout/dropout/GreaterEqual/y?
4sequential/functional_1/dropout/dropout/GreaterEqualGreaterEqualMsequential/functional_1/dropout/dropout/random_uniform/RandomUniform:output:0?sequential/functional_1/dropout/dropout/GreaterEqual/y:output:0*
T0*(
_output_shapes
:??????????
26
4sequential/functional_1/dropout/dropout/GreaterEqual?
,sequential/functional_1/dropout/dropout/CastCast8sequential/functional_1/dropout/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*(
_output_shapes
:??????????
2.
,sequential/functional_1/dropout/dropout/Cast?
-sequential/functional_1/dropout/dropout/Mul_1Mul/sequential/functional_1/dropout/dropout/Mul:z:00sequential/functional_1/dropout/dropout/Cast:y:0*
T0*(
_output_shapes
:??????????
2/
-sequential/functional_1/dropout/dropout/Mul_1?
4sequential/functional_1/layer1/MatMul/ReadVariableOpReadVariableOp=sequential_functional_1_layer1_matmul_readvariableop_resource* 
_output_shapes
:
?
?*
dtype026
4sequential/functional_1/layer1/MatMul/ReadVariableOp?
%sequential/functional_1/layer1/MatMulMatMul1sequential/functional_1/dropout/dropout/Mul_1:z:0<sequential/functional_1/layer1/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2'
%sequential/functional_1/layer1/MatMul?
5sequential/functional_1/layer1/BiasAdd/ReadVariableOpReadVariableOp>sequential_functional_1_layer1_biasadd_readvariableop_resource*
_output_shapes	
:?*
dtype027
5sequential/functional_1/layer1/BiasAdd/ReadVariableOp?
&sequential/functional_1/layer1/BiasAddBiasAdd/sequential/functional_1/layer1/MatMul:product:0=sequential/functional_1/layer1/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2(
&sequential/functional_1/layer1/BiasAdd?
"sequential/functional_1/layer1/EluElu/sequential/functional_1/layer1/BiasAdd:output:0*
T0*(
_output_shapes
:??????????2$
"sequential/functional_1/layer1/Elu?
&sequential/dense/MatMul/ReadVariableOpReadVariableOp/sequential_dense_matmul_readvariableop_resource*
_output_shapes
:	?d*
dtype02(
&sequential/dense/MatMul/ReadVariableOp?
sequential/dense/MatMulMatMul0sequential/functional_1/layer1/Elu:activations:0.sequential/dense/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????d2
sequential/dense/MatMul?
'sequential/dense/BiasAdd/ReadVariableOpReadVariableOp0sequential_dense_biasadd_readvariableop_resource*
_output_shapes
:d*
dtype02)
'sequential/dense/BiasAdd/ReadVariableOp?
sequential/dense/BiasAddBiasAdd!sequential/dense/MatMul:product:0/sequential/dense/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????d2
sequential/dense/BiasAdd?
sequential/dense/EluElu!sequential/dense/BiasAdd:output:0*
T0*'
_output_shapes
:?????????d2
sequential/dense/Elu?
(sequential/dense_1/MatMul/ReadVariableOpReadVariableOp1sequential_dense_1_matmul_readvariableop_resource*
_output_shapes

:d
*
dtype02*
(sequential/dense_1/MatMul/ReadVariableOp?
sequential/dense_1/MatMulMatMul"sequential/dense/Elu:activations:00sequential/dense_1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????
2
sequential/dense_1/MatMul?
)sequential/dense_1/BiasAdd/ReadVariableOpReadVariableOp2sequential_dense_1_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02+
)sequential/dense_1/BiasAdd/ReadVariableOp?
sequential/dense_1/BiasAddBiasAdd#sequential/dense_1/MatMul:product:01sequential/dense_1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????
2
sequential/dense_1/BiasAdd?
sequential/dense_1/EluElu#sequential/dense_1/BiasAdd:output:0*
T0*'
_output_shapes
:?????????
2
sequential/dense_1/Elu?
6sequential/functional_1/conv3d/Conv3D_1/ReadVariableOpReadVariableOp=sequential_functional_1_conv3d_conv3d_readvariableop_resource*+
_output_shapes
:?*
dtype028
6sequential/functional_1/conv3d/Conv3D_1/ReadVariableOp?
'sequential/functional_1/conv3d/Conv3D_1Conv3Dinputs_1>sequential/functional_1/conv3d/Conv3D_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2)
'sequential/functional_1/conv3d/Conv3D_1?
7sequential/functional_1/conv3d/BiasAdd_1/ReadVariableOpReadVariableOp>sequential_functional_1_conv3d_biasadd_readvariableop_resource*
_output_shapes
:*
dtype029
7sequential/functional_1/conv3d/BiasAdd_1/ReadVariableOp?
(sequential/functional_1/conv3d/BiasAdd_1BiasAdd0sequential/functional_1/conv3d/Conv3D_1:output:0?sequential/functional_1/conv3d/BiasAdd_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2*
(sequential/functional_1/conv3d/BiasAdd_1?
Lsequential/functional_1/batch_normalization/moments_1/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2N
Lsequential/functional_1/batch_normalization/moments_1/mean/reduction_indices?
:sequential/functional_1/batch_normalization/moments_1/meanMean1sequential/functional_1/conv3d/BiasAdd_1:output:0Usequential/functional_1/batch_normalization/moments_1/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2<
:sequential/functional_1/batch_normalization/moments_1/mean?
Bsequential/functional_1/batch_normalization/moments_1/StopGradientStopGradientCsequential/functional_1/batch_normalization/moments_1/mean:output:0*
T0**
_output_shapes
:2D
Bsequential/functional_1/batch_normalization/moments_1/StopGradient?
Gsequential/functional_1/batch_normalization/moments_1/SquaredDifferenceSquaredDifference1sequential/functional_1/conv3d/BiasAdd_1:output:0Ksequential/functional_1/batch_normalization/moments_1/StopGradient:output:0*
T0*3
_output_shapes!
:?????????2I
Gsequential/functional_1/batch_normalization/moments_1/SquaredDifference?
Psequential/functional_1/batch_normalization/moments_1/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2R
Psequential/functional_1/batch_normalization/moments_1/variance/reduction_indices?
>sequential/functional_1/batch_normalization/moments_1/varianceMeanKsequential/functional_1/batch_normalization/moments_1/SquaredDifference:z:0Ysequential/functional_1/batch_normalization/moments_1/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2@
>sequential/functional_1/batch_normalization/moments_1/variance?
=sequential/functional_1/batch_normalization/moments_1/SqueezeSqueezeCsequential/functional_1/batch_normalization/moments_1/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2?
=sequential/functional_1/batch_normalization/moments_1/Squeeze?
?sequential/functional_1/batch_normalization/moments_1/Squeeze_1SqueezeGsequential/functional_1/batch_normalization/moments_1/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2A
?sequential/functional_1/batch_normalization/moments_1/Squeeze_1?
Csequential/functional_1/batch_normalization/AssignMovingAvg_2/decayConst*V
_classL
JHloc:@sequential/functional_1/batch_normalization/AssignMovingAvg/1225167*
_output_shapes
: *
dtype0*
valueB
 *
?#<2E
Csequential/functional_1/batch_normalization/AssignMovingAvg_2/decay?
Lsequential/functional_1/batch_normalization/AssignMovingAvg_2/ReadVariableOpReadVariableOpCsequential_functional_1_batch_normalization_assignmovingavg_1225167P^sequential/functional_1/batch_normalization/AssignMovingAvg/AssignSubVariableOp*
_output_shapes
:*
dtype02N
Lsequential/functional_1/batch_normalization/AssignMovingAvg_2/ReadVariableOp?
Asequential/functional_1/batch_normalization/AssignMovingAvg_2/subSubTsequential/functional_1/batch_normalization/AssignMovingAvg_2/ReadVariableOp:value:0Fsequential/functional_1/batch_normalization/moments_1/Squeeze:output:0*
T0*V
_classL
JHloc:@sequential/functional_1/batch_normalization/AssignMovingAvg/1225167*
_output_shapes
:2C
Asequential/functional_1/batch_normalization/AssignMovingAvg_2/sub?
Asequential/functional_1/batch_normalization/AssignMovingAvg_2/mulMulEsequential/functional_1/batch_normalization/AssignMovingAvg_2/sub:z:0Lsequential/functional_1/batch_normalization/AssignMovingAvg_2/decay:output:0*
T0*V
_classL
JHloc:@sequential/functional_1/batch_normalization/AssignMovingAvg/1225167*
_output_shapes
:2C
Asequential/functional_1/batch_normalization/AssignMovingAvg_2/mul?
Qsequential/functional_1/batch_normalization/AssignMovingAvg_2/AssignSubVariableOpAssignSubVariableOpCsequential_functional_1_batch_normalization_assignmovingavg_1225167Esequential/functional_1/batch_normalization/AssignMovingAvg_2/mul:z:0P^sequential/functional_1/batch_normalization/AssignMovingAvg/AssignSubVariableOpM^sequential/functional_1/batch_normalization/AssignMovingAvg_2/ReadVariableOp*V
_classL
JHloc:@sequential/functional_1/batch_normalization/AssignMovingAvg/1225167*
_output_shapes
 *
dtype02S
Qsequential/functional_1/batch_normalization/AssignMovingAvg_2/AssignSubVariableOp?
Csequential/functional_1/batch_normalization/AssignMovingAvg_3/decayConst*X
_classN
LJloc:@sequential/functional_1/batch_normalization/AssignMovingAvg_1/1225173*
_output_shapes
: *
dtype0*
valueB
 *
?#<2E
Csequential/functional_1/batch_normalization/AssignMovingAvg_3/decay?
Lsequential/functional_1/batch_normalization/AssignMovingAvg_3/ReadVariableOpReadVariableOpEsequential_functional_1_batch_normalization_assignmovingavg_1_1225173R^sequential/functional_1/batch_normalization/AssignMovingAvg_1/AssignSubVariableOp*
_output_shapes
:*
dtype02N
Lsequential/functional_1/batch_normalization/AssignMovingAvg_3/ReadVariableOp?
Asequential/functional_1/batch_normalization/AssignMovingAvg_3/subSubTsequential/functional_1/batch_normalization/AssignMovingAvg_3/ReadVariableOp:value:0Hsequential/functional_1/batch_normalization/moments_1/Squeeze_1:output:0*
T0*X
_classN
LJloc:@sequential/functional_1/batch_normalization/AssignMovingAvg_1/1225173*
_output_shapes
:2C
Asequential/functional_1/batch_normalization/AssignMovingAvg_3/sub?
Asequential/functional_1/batch_normalization/AssignMovingAvg_3/mulMulEsequential/functional_1/batch_normalization/AssignMovingAvg_3/sub:z:0Lsequential/functional_1/batch_normalization/AssignMovingAvg_3/decay:output:0*
T0*X
_classN
LJloc:@sequential/functional_1/batch_normalization/AssignMovingAvg_1/1225173*
_output_shapes
:2C
Asequential/functional_1/batch_normalization/AssignMovingAvg_3/mul?
Qsequential/functional_1/batch_normalization/AssignMovingAvg_3/AssignSubVariableOpAssignSubVariableOpEsequential_functional_1_batch_normalization_assignmovingavg_1_1225173Esequential/functional_1/batch_normalization/AssignMovingAvg_3/mul:z:0R^sequential/functional_1/batch_normalization/AssignMovingAvg_1/AssignSubVariableOpM^sequential/functional_1/batch_normalization/AssignMovingAvg_3/ReadVariableOp*X
_classN
LJloc:@sequential/functional_1/batch_normalization/AssignMovingAvg_1/1225173*
_output_shapes
 *
dtype02S
Qsequential/functional_1/batch_normalization/AssignMovingAvg_3/AssignSubVariableOp?
=sequential/functional_1/batch_normalization/batchnorm_1/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2?
=sequential/functional_1/batch_normalization/batchnorm_1/add/y?
;sequential/functional_1/batch_normalization/batchnorm_1/addAddV2Hsequential/functional_1/batch_normalization/moments_1/Squeeze_1:output:0Fsequential/functional_1/batch_normalization/batchnorm_1/add/y:output:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization/batchnorm_1/add?
=sequential/functional_1/batch_normalization/batchnorm_1/RsqrtRsqrt?sequential/functional_1/batch_normalization/batchnorm_1/add:z:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization/batchnorm_1/Rsqrt?
Jsequential/functional_1/batch_normalization/batchnorm_1/mul/ReadVariableOpReadVariableOpQsequential_functional_1_batch_normalization_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02L
Jsequential/functional_1/batch_normalization/batchnorm_1/mul/ReadVariableOp?
;sequential/functional_1/batch_normalization/batchnorm_1/mulMulAsequential/functional_1/batch_normalization/batchnorm_1/Rsqrt:y:0Rsequential/functional_1/batch_normalization/batchnorm_1/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization/batchnorm_1/mul?
=sequential/functional_1/batch_normalization/batchnorm_1/mul_1Mul1sequential/functional_1/conv3d/BiasAdd_1:output:0?sequential/functional_1/batch_normalization/batchnorm_1/mul:z:0*
T0*3
_output_shapes!
:?????????2?
=sequential/functional_1/batch_normalization/batchnorm_1/mul_1?
=sequential/functional_1/batch_normalization/batchnorm_1/mul_2MulFsequential/functional_1/batch_normalization/moments_1/Squeeze:output:0?sequential/functional_1/batch_normalization/batchnorm_1/mul:z:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization/batchnorm_1/mul_2?
Fsequential/functional_1/batch_normalization/batchnorm_1/ReadVariableOpReadVariableOpMsequential_functional_1_batch_normalization_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02H
Fsequential/functional_1/batch_normalization/batchnorm_1/ReadVariableOp?
;sequential/functional_1/batch_normalization/batchnorm_1/subSubNsequential/functional_1/batch_normalization/batchnorm_1/ReadVariableOp:value:0Asequential/functional_1/batch_normalization/batchnorm_1/mul_2:z:0*
T0*
_output_shapes
:2=
;sequential/functional_1/batch_normalization/batchnorm_1/sub?
=sequential/functional_1/batch_normalization/batchnorm_1/add_1AddV2Asequential/functional_1/batch_normalization/batchnorm_1/mul_1:z:0?sequential/functional_1/batch_normalization/batchnorm_1/sub:z:0*
T0*3
_output_shapes!
:?????????2?
=sequential/functional_1/batch_normalization/batchnorm_1/add_1?
8sequential/functional_1/conv3d_1/Conv3D_1/ReadVariableOpReadVariableOp?sequential_functional_1_conv3d_1_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02:
8sequential/functional_1/conv3d_1/Conv3D_1/ReadVariableOp?
)sequential/functional_1/conv3d_1/Conv3D_1Conv3DAsequential/functional_1/batch_normalization/batchnorm_1/add_1:z:0@sequential/functional_1/conv3d_1/Conv3D_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2+
)sequential/functional_1/conv3d_1/Conv3D_1?
9sequential/functional_1/conv3d_1/BiasAdd_1/ReadVariableOpReadVariableOp@sequential_functional_1_conv3d_1_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02;
9sequential/functional_1/conv3d_1/BiasAdd_1/ReadVariableOp?
*sequential/functional_1/conv3d_1/BiasAdd_1BiasAdd2sequential/functional_1/conv3d_1/Conv3D_1:output:0Asequential/functional_1/conv3d_1/BiasAdd_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2,
*sequential/functional_1/conv3d_1/BiasAdd_1?
&sequential/functional_1/conv3d_1/Elu_1Elu3sequential/functional_1/conv3d_1/BiasAdd_1:output:0*
T0*3
_output_shapes!
:?????????2(
&sequential/functional_1/conv3d_1/Elu_1?
Nsequential/functional_1/batch_normalization_1/moments_1/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2P
Nsequential/functional_1/batch_normalization_1/moments_1/mean/reduction_indices?
<sequential/functional_1/batch_normalization_1/moments_1/meanMean4sequential/functional_1/conv3d_1/Elu_1:activations:0Wsequential/functional_1/batch_normalization_1/moments_1/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2>
<sequential/functional_1/batch_normalization_1/moments_1/mean?
Dsequential/functional_1/batch_normalization_1/moments_1/StopGradientStopGradientEsequential/functional_1/batch_normalization_1/moments_1/mean:output:0*
T0**
_output_shapes
:2F
Dsequential/functional_1/batch_normalization_1/moments_1/StopGradient?
Isequential/functional_1/batch_normalization_1/moments_1/SquaredDifferenceSquaredDifference4sequential/functional_1/conv3d_1/Elu_1:activations:0Msequential/functional_1/batch_normalization_1/moments_1/StopGradient:output:0*
T0*3
_output_shapes!
:?????????2K
Isequential/functional_1/batch_normalization_1/moments_1/SquaredDifference?
Rsequential/functional_1/batch_normalization_1/moments_1/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2T
Rsequential/functional_1/batch_normalization_1/moments_1/variance/reduction_indices?
@sequential/functional_1/batch_normalization_1/moments_1/varianceMeanMsequential/functional_1/batch_normalization_1/moments_1/SquaredDifference:z:0[sequential/functional_1/batch_normalization_1/moments_1/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2B
@sequential/functional_1/batch_normalization_1/moments_1/variance?
?sequential/functional_1/batch_normalization_1/moments_1/SqueezeSqueezeEsequential/functional_1/batch_normalization_1/moments_1/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2A
?sequential/functional_1/batch_normalization_1/moments_1/Squeeze?
Asequential/functional_1/batch_normalization_1/moments_1/Squeeze_1SqueezeIsequential/functional_1/batch_normalization_1/moments_1/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2C
Asequential/functional_1/batch_normalization_1/moments_1/Squeeze_1?
Esequential/functional_1/batch_normalization_1/AssignMovingAvg_2/decayConst*X
_classN
LJloc:@sequential/functional_1/batch_normalization_1/AssignMovingAvg/1225206*
_output_shapes
: *
dtype0*
valueB
 *
?#<2G
Esequential/functional_1/batch_normalization_1/AssignMovingAvg_2/decay?
Nsequential/functional_1/batch_normalization_1/AssignMovingAvg_2/ReadVariableOpReadVariableOpEsequential_functional_1_batch_normalization_1_assignmovingavg_1225206R^sequential/functional_1/batch_normalization_1/AssignMovingAvg/AssignSubVariableOp*
_output_shapes
:*
dtype02P
Nsequential/functional_1/batch_normalization_1/AssignMovingAvg_2/ReadVariableOp?
Csequential/functional_1/batch_normalization_1/AssignMovingAvg_2/subSubVsequential/functional_1/batch_normalization_1/AssignMovingAvg_2/ReadVariableOp:value:0Hsequential/functional_1/batch_normalization_1/moments_1/Squeeze:output:0*
T0*X
_classN
LJloc:@sequential/functional_1/batch_normalization_1/AssignMovingAvg/1225206*
_output_shapes
:2E
Csequential/functional_1/batch_normalization_1/AssignMovingAvg_2/sub?
Csequential/functional_1/batch_normalization_1/AssignMovingAvg_2/mulMulGsequential/functional_1/batch_normalization_1/AssignMovingAvg_2/sub:z:0Nsequential/functional_1/batch_normalization_1/AssignMovingAvg_2/decay:output:0*
T0*X
_classN
LJloc:@sequential/functional_1/batch_normalization_1/AssignMovingAvg/1225206*
_output_shapes
:2E
Csequential/functional_1/batch_normalization_1/AssignMovingAvg_2/mul?
Ssequential/functional_1/batch_normalization_1/AssignMovingAvg_2/AssignSubVariableOpAssignSubVariableOpEsequential_functional_1_batch_normalization_1_assignmovingavg_1225206Gsequential/functional_1/batch_normalization_1/AssignMovingAvg_2/mul:z:0R^sequential/functional_1/batch_normalization_1/AssignMovingAvg/AssignSubVariableOpO^sequential/functional_1/batch_normalization_1/AssignMovingAvg_2/ReadVariableOp*X
_classN
LJloc:@sequential/functional_1/batch_normalization_1/AssignMovingAvg/1225206*
_output_shapes
 *
dtype02U
Ssequential/functional_1/batch_normalization_1/AssignMovingAvg_2/AssignSubVariableOp?
Esequential/functional_1/batch_normalization_1/AssignMovingAvg_3/decayConst*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_1/AssignMovingAvg_1/1225212*
_output_shapes
: *
dtype0*
valueB
 *
?#<2G
Esequential/functional_1/batch_normalization_1/AssignMovingAvg_3/decay?
Nsequential/functional_1/batch_normalization_1/AssignMovingAvg_3/ReadVariableOpReadVariableOpGsequential_functional_1_batch_normalization_1_assignmovingavg_1_1225212T^sequential/functional_1/batch_normalization_1/AssignMovingAvg_1/AssignSubVariableOp*
_output_shapes
:*
dtype02P
Nsequential/functional_1/batch_normalization_1/AssignMovingAvg_3/ReadVariableOp?
Csequential/functional_1/batch_normalization_1/AssignMovingAvg_3/subSubVsequential/functional_1/batch_normalization_1/AssignMovingAvg_3/ReadVariableOp:value:0Jsequential/functional_1/batch_normalization_1/moments_1/Squeeze_1:output:0*
T0*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_1/AssignMovingAvg_1/1225212*
_output_shapes
:2E
Csequential/functional_1/batch_normalization_1/AssignMovingAvg_3/sub?
Csequential/functional_1/batch_normalization_1/AssignMovingAvg_3/mulMulGsequential/functional_1/batch_normalization_1/AssignMovingAvg_3/sub:z:0Nsequential/functional_1/batch_normalization_1/AssignMovingAvg_3/decay:output:0*
T0*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_1/AssignMovingAvg_1/1225212*
_output_shapes
:2E
Csequential/functional_1/batch_normalization_1/AssignMovingAvg_3/mul?
Ssequential/functional_1/batch_normalization_1/AssignMovingAvg_3/AssignSubVariableOpAssignSubVariableOpGsequential_functional_1_batch_normalization_1_assignmovingavg_1_1225212Gsequential/functional_1/batch_normalization_1/AssignMovingAvg_3/mul:z:0T^sequential/functional_1/batch_normalization_1/AssignMovingAvg_1/AssignSubVariableOpO^sequential/functional_1/batch_normalization_1/AssignMovingAvg_3/ReadVariableOp*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_1/AssignMovingAvg_1/1225212*
_output_shapes
 *
dtype02U
Ssequential/functional_1/batch_normalization_1/AssignMovingAvg_3/AssignSubVariableOp?
?sequential/functional_1/batch_normalization_1/batchnorm_1/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2A
?sequential/functional_1/batch_normalization_1/batchnorm_1/add/y?
=sequential/functional_1/batch_normalization_1/batchnorm_1/addAddV2Jsequential/functional_1/batch_normalization_1/moments_1/Squeeze_1:output:0Hsequential/functional_1/batch_normalization_1/batchnorm_1/add/y:output:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_1/batchnorm_1/add?
?sequential/functional_1/batch_normalization_1/batchnorm_1/RsqrtRsqrtAsequential/functional_1/batch_normalization_1/batchnorm_1/add:z:0*
T0*
_output_shapes
:2A
?sequential/functional_1/batch_normalization_1/batchnorm_1/Rsqrt?
Lsequential/functional_1/batch_normalization_1/batchnorm_1/mul/ReadVariableOpReadVariableOpSsequential_functional_1_batch_normalization_1_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02N
Lsequential/functional_1/batch_normalization_1/batchnorm_1/mul/ReadVariableOp?
=sequential/functional_1/batch_normalization_1/batchnorm_1/mulMulCsequential/functional_1/batch_normalization_1/batchnorm_1/Rsqrt:y:0Tsequential/functional_1/batch_normalization_1/batchnorm_1/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_1/batchnorm_1/mul?
?sequential/functional_1/batch_normalization_1/batchnorm_1/mul_1Mul4sequential/functional_1/conv3d_1/Elu_1:activations:0Asequential/functional_1/batch_normalization_1/batchnorm_1/mul:z:0*
T0*3
_output_shapes!
:?????????2A
?sequential/functional_1/batch_normalization_1/batchnorm_1/mul_1?
?sequential/functional_1/batch_normalization_1/batchnorm_1/mul_2MulHsequential/functional_1/batch_normalization_1/moments_1/Squeeze:output:0Asequential/functional_1/batch_normalization_1/batchnorm_1/mul:z:0*
T0*
_output_shapes
:2A
?sequential/functional_1/batch_normalization_1/batchnorm_1/mul_2?
Hsequential/functional_1/batch_normalization_1/batchnorm_1/ReadVariableOpReadVariableOpOsequential_functional_1_batch_normalization_1_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02J
Hsequential/functional_1/batch_normalization_1/batchnorm_1/ReadVariableOp?
=sequential/functional_1/batch_normalization_1/batchnorm_1/subSubPsequential/functional_1/batch_normalization_1/batchnorm_1/ReadVariableOp:value:0Csequential/functional_1/batch_normalization_1/batchnorm_1/mul_2:z:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_1/batchnorm_1/sub?
?sequential/functional_1/batch_normalization_1/batchnorm_1/add_1AddV2Csequential/functional_1/batch_normalization_1/batchnorm_1/mul_1:z:0Asequential/functional_1/batch_normalization_1/batchnorm_1/sub:z:0*
T0*3
_output_shapes!
:?????????2A
?sequential/functional_1/batch_normalization_1/batchnorm_1/add_1?
8sequential/functional_1/conv3d_2/Conv3D_1/ReadVariableOpReadVariableOp?sequential_functional_1_conv3d_2_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02:
8sequential/functional_1/conv3d_2/Conv3D_1/ReadVariableOp?
)sequential/functional_1/conv3d_2/Conv3D_1Conv3DCsequential/functional_1/batch_normalization_1/batchnorm_1/add_1:z:0@sequential/functional_1/conv3d_2/Conv3D_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2+
)sequential/functional_1/conv3d_2/Conv3D_1?
9sequential/functional_1/conv3d_2/BiasAdd_1/ReadVariableOpReadVariableOp@sequential_functional_1_conv3d_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02;
9sequential/functional_1/conv3d_2/BiasAdd_1/ReadVariableOp?
*sequential/functional_1/conv3d_2/BiasAdd_1BiasAdd2sequential/functional_1/conv3d_2/Conv3D_1:output:0Asequential/functional_1/conv3d_2/BiasAdd_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2,
*sequential/functional_1/conv3d_2/BiasAdd_1?
&sequential/functional_1/conv3d_2/Elu_1Elu3sequential/functional_1/conv3d_2/BiasAdd_1:output:0*
T0*3
_output_shapes!
:?????????2(
&sequential/functional_1/conv3d_2/Elu_1?
Nsequential/functional_1/batch_normalization_2/moments_1/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2P
Nsequential/functional_1/batch_normalization_2/moments_1/mean/reduction_indices?
<sequential/functional_1/batch_normalization_2/moments_1/meanMean4sequential/functional_1/conv3d_2/Elu_1:activations:0Wsequential/functional_1/batch_normalization_2/moments_1/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2>
<sequential/functional_1/batch_normalization_2/moments_1/mean?
Dsequential/functional_1/batch_normalization_2/moments_1/StopGradientStopGradientEsequential/functional_1/batch_normalization_2/moments_1/mean:output:0*
T0**
_output_shapes
:2F
Dsequential/functional_1/batch_normalization_2/moments_1/StopGradient?
Isequential/functional_1/batch_normalization_2/moments_1/SquaredDifferenceSquaredDifference4sequential/functional_1/conv3d_2/Elu_1:activations:0Msequential/functional_1/batch_normalization_2/moments_1/StopGradient:output:0*
T0*3
_output_shapes!
:?????????2K
Isequential/functional_1/batch_normalization_2/moments_1/SquaredDifference?
Rsequential/functional_1/batch_normalization_2/moments_1/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2T
Rsequential/functional_1/batch_normalization_2/moments_1/variance/reduction_indices?
@sequential/functional_1/batch_normalization_2/moments_1/varianceMeanMsequential/functional_1/batch_normalization_2/moments_1/SquaredDifference:z:0[sequential/functional_1/batch_normalization_2/moments_1/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2B
@sequential/functional_1/batch_normalization_2/moments_1/variance?
?sequential/functional_1/batch_normalization_2/moments_1/SqueezeSqueezeEsequential/functional_1/batch_normalization_2/moments_1/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2A
?sequential/functional_1/batch_normalization_2/moments_1/Squeeze?
Asequential/functional_1/batch_normalization_2/moments_1/Squeeze_1SqueezeIsequential/functional_1/batch_normalization_2/moments_1/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2C
Asequential/functional_1/batch_normalization_2/moments_1/Squeeze_1?
Esequential/functional_1/batch_normalization_2/AssignMovingAvg_2/decayConst*X
_classN
LJloc:@sequential/functional_1/batch_normalization_2/AssignMovingAvg/1225245*
_output_shapes
: *
dtype0*
valueB
 *
?#<2G
Esequential/functional_1/batch_normalization_2/AssignMovingAvg_2/decay?
Nsequential/functional_1/batch_normalization_2/AssignMovingAvg_2/ReadVariableOpReadVariableOpEsequential_functional_1_batch_normalization_2_assignmovingavg_1225245R^sequential/functional_1/batch_normalization_2/AssignMovingAvg/AssignSubVariableOp*
_output_shapes
:*
dtype02P
Nsequential/functional_1/batch_normalization_2/AssignMovingAvg_2/ReadVariableOp?
Csequential/functional_1/batch_normalization_2/AssignMovingAvg_2/subSubVsequential/functional_1/batch_normalization_2/AssignMovingAvg_2/ReadVariableOp:value:0Hsequential/functional_1/batch_normalization_2/moments_1/Squeeze:output:0*
T0*X
_classN
LJloc:@sequential/functional_1/batch_normalization_2/AssignMovingAvg/1225245*
_output_shapes
:2E
Csequential/functional_1/batch_normalization_2/AssignMovingAvg_2/sub?
Csequential/functional_1/batch_normalization_2/AssignMovingAvg_2/mulMulGsequential/functional_1/batch_normalization_2/AssignMovingAvg_2/sub:z:0Nsequential/functional_1/batch_normalization_2/AssignMovingAvg_2/decay:output:0*
T0*X
_classN
LJloc:@sequential/functional_1/batch_normalization_2/AssignMovingAvg/1225245*
_output_shapes
:2E
Csequential/functional_1/batch_normalization_2/AssignMovingAvg_2/mul?
Ssequential/functional_1/batch_normalization_2/AssignMovingAvg_2/AssignSubVariableOpAssignSubVariableOpEsequential_functional_1_batch_normalization_2_assignmovingavg_1225245Gsequential/functional_1/batch_normalization_2/AssignMovingAvg_2/mul:z:0R^sequential/functional_1/batch_normalization_2/AssignMovingAvg/AssignSubVariableOpO^sequential/functional_1/batch_normalization_2/AssignMovingAvg_2/ReadVariableOp*X
_classN
LJloc:@sequential/functional_1/batch_normalization_2/AssignMovingAvg/1225245*
_output_shapes
 *
dtype02U
Ssequential/functional_1/batch_normalization_2/AssignMovingAvg_2/AssignSubVariableOp?
Esequential/functional_1/batch_normalization_2/AssignMovingAvg_3/decayConst*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_2/AssignMovingAvg_1/1225251*
_output_shapes
: *
dtype0*
valueB
 *
?#<2G
Esequential/functional_1/batch_normalization_2/AssignMovingAvg_3/decay?
Nsequential/functional_1/batch_normalization_2/AssignMovingAvg_3/ReadVariableOpReadVariableOpGsequential_functional_1_batch_normalization_2_assignmovingavg_1_1225251T^sequential/functional_1/batch_normalization_2/AssignMovingAvg_1/AssignSubVariableOp*
_output_shapes
:*
dtype02P
Nsequential/functional_1/batch_normalization_2/AssignMovingAvg_3/ReadVariableOp?
Csequential/functional_1/batch_normalization_2/AssignMovingAvg_3/subSubVsequential/functional_1/batch_normalization_2/AssignMovingAvg_3/ReadVariableOp:value:0Jsequential/functional_1/batch_normalization_2/moments_1/Squeeze_1:output:0*
T0*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_2/AssignMovingAvg_1/1225251*
_output_shapes
:2E
Csequential/functional_1/batch_normalization_2/AssignMovingAvg_3/sub?
Csequential/functional_1/batch_normalization_2/AssignMovingAvg_3/mulMulGsequential/functional_1/batch_normalization_2/AssignMovingAvg_3/sub:z:0Nsequential/functional_1/batch_normalization_2/AssignMovingAvg_3/decay:output:0*
T0*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_2/AssignMovingAvg_1/1225251*
_output_shapes
:2E
Csequential/functional_1/batch_normalization_2/AssignMovingAvg_3/mul?
Ssequential/functional_1/batch_normalization_2/AssignMovingAvg_3/AssignSubVariableOpAssignSubVariableOpGsequential_functional_1_batch_normalization_2_assignmovingavg_1_1225251Gsequential/functional_1/batch_normalization_2/AssignMovingAvg_3/mul:z:0T^sequential/functional_1/batch_normalization_2/AssignMovingAvg_1/AssignSubVariableOpO^sequential/functional_1/batch_normalization_2/AssignMovingAvg_3/ReadVariableOp*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_2/AssignMovingAvg_1/1225251*
_output_shapes
 *
dtype02U
Ssequential/functional_1/batch_normalization_2/AssignMovingAvg_3/AssignSubVariableOp?
?sequential/functional_1/batch_normalization_2/batchnorm_1/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2A
?sequential/functional_1/batch_normalization_2/batchnorm_1/add/y?
=sequential/functional_1/batch_normalization_2/batchnorm_1/addAddV2Jsequential/functional_1/batch_normalization_2/moments_1/Squeeze_1:output:0Hsequential/functional_1/batch_normalization_2/batchnorm_1/add/y:output:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_2/batchnorm_1/add?
?sequential/functional_1/batch_normalization_2/batchnorm_1/RsqrtRsqrtAsequential/functional_1/batch_normalization_2/batchnorm_1/add:z:0*
T0*
_output_shapes
:2A
?sequential/functional_1/batch_normalization_2/batchnorm_1/Rsqrt?
Lsequential/functional_1/batch_normalization_2/batchnorm_1/mul/ReadVariableOpReadVariableOpSsequential_functional_1_batch_normalization_2_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02N
Lsequential/functional_1/batch_normalization_2/batchnorm_1/mul/ReadVariableOp?
=sequential/functional_1/batch_normalization_2/batchnorm_1/mulMulCsequential/functional_1/batch_normalization_2/batchnorm_1/Rsqrt:y:0Tsequential/functional_1/batch_normalization_2/batchnorm_1/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_2/batchnorm_1/mul?
?sequential/functional_1/batch_normalization_2/batchnorm_1/mul_1Mul4sequential/functional_1/conv3d_2/Elu_1:activations:0Asequential/functional_1/batch_normalization_2/batchnorm_1/mul:z:0*
T0*3
_output_shapes!
:?????????2A
?sequential/functional_1/batch_normalization_2/batchnorm_1/mul_1?
?sequential/functional_1/batch_normalization_2/batchnorm_1/mul_2MulHsequential/functional_1/batch_normalization_2/moments_1/Squeeze:output:0Asequential/functional_1/batch_normalization_2/batchnorm_1/mul:z:0*
T0*
_output_shapes
:2A
?sequential/functional_1/batch_normalization_2/batchnorm_1/mul_2?
Hsequential/functional_1/batch_normalization_2/batchnorm_1/ReadVariableOpReadVariableOpOsequential_functional_1_batch_normalization_2_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02J
Hsequential/functional_1/batch_normalization_2/batchnorm_1/ReadVariableOp?
=sequential/functional_1/batch_normalization_2/batchnorm_1/subSubPsequential/functional_1/batch_normalization_2/batchnorm_1/ReadVariableOp:value:0Csequential/functional_1/batch_normalization_2/batchnorm_1/mul_2:z:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_2/batchnorm_1/sub?
?sequential/functional_1/batch_normalization_2/batchnorm_1/add_1AddV2Csequential/functional_1/batch_normalization_2/batchnorm_1/mul_1:z:0Asequential/functional_1/batch_normalization_2/batchnorm_1/sub:z:0*
T0*3
_output_shapes!
:?????????2A
?sequential/functional_1/batch_normalization_2/batchnorm_1/add_1?
8sequential/functional_1/conv3d_3/Conv3D_1/ReadVariableOpReadVariableOp?sequential_functional_1_conv3d_3_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02:
8sequential/functional_1/conv3d_3/Conv3D_1/ReadVariableOp?
)sequential/functional_1/conv3d_3/Conv3D_1Conv3DCsequential/functional_1/batch_normalization_2/batchnorm_1/add_1:z:0@sequential/functional_1/conv3d_3/Conv3D_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2+
)sequential/functional_1/conv3d_3/Conv3D_1?
9sequential/functional_1/conv3d_3/BiasAdd_1/ReadVariableOpReadVariableOp@sequential_functional_1_conv3d_3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02;
9sequential/functional_1/conv3d_3/BiasAdd_1/ReadVariableOp?
*sequential/functional_1/conv3d_3/BiasAdd_1BiasAdd2sequential/functional_1/conv3d_3/Conv3D_1:output:0Asequential/functional_1/conv3d_3/BiasAdd_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2,
*sequential/functional_1/conv3d_3/BiasAdd_1?
&sequential/functional_1/conv3d_3/Elu_1Elu3sequential/functional_1/conv3d_3/BiasAdd_1:output:0*
T0*3
_output_shapes!
:?????????2(
&sequential/functional_1/conv3d_3/Elu_1?
Nsequential/functional_1/batch_normalization_3/moments_1/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2P
Nsequential/functional_1/batch_normalization_3/moments_1/mean/reduction_indices?
<sequential/functional_1/batch_normalization_3/moments_1/meanMean4sequential/functional_1/conv3d_3/Elu_1:activations:0Wsequential/functional_1/batch_normalization_3/moments_1/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2>
<sequential/functional_1/batch_normalization_3/moments_1/mean?
Dsequential/functional_1/batch_normalization_3/moments_1/StopGradientStopGradientEsequential/functional_1/batch_normalization_3/moments_1/mean:output:0*
T0**
_output_shapes
:2F
Dsequential/functional_1/batch_normalization_3/moments_1/StopGradient?
Isequential/functional_1/batch_normalization_3/moments_1/SquaredDifferenceSquaredDifference4sequential/functional_1/conv3d_3/Elu_1:activations:0Msequential/functional_1/batch_normalization_3/moments_1/StopGradient:output:0*
T0*3
_output_shapes!
:?????????2K
Isequential/functional_1/batch_normalization_3/moments_1/SquaredDifference?
Rsequential/functional_1/batch_normalization_3/moments_1/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2T
Rsequential/functional_1/batch_normalization_3/moments_1/variance/reduction_indices?
@sequential/functional_1/batch_normalization_3/moments_1/varianceMeanMsequential/functional_1/batch_normalization_3/moments_1/SquaredDifference:z:0[sequential/functional_1/batch_normalization_3/moments_1/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2B
@sequential/functional_1/batch_normalization_3/moments_1/variance?
?sequential/functional_1/batch_normalization_3/moments_1/SqueezeSqueezeEsequential/functional_1/batch_normalization_3/moments_1/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2A
?sequential/functional_1/batch_normalization_3/moments_1/Squeeze?
Asequential/functional_1/batch_normalization_3/moments_1/Squeeze_1SqueezeIsequential/functional_1/batch_normalization_3/moments_1/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2C
Asequential/functional_1/batch_normalization_3/moments_1/Squeeze_1?
Esequential/functional_1/batch_normalization_3/AssignMovingAvg_2/decayConst*X
_classN
LJloc:@sequential/functional_1/batch_normalization_3/AssignMovingAvg/1225284*
_output_shapes
: *
dtype0*
valueB
 *
?#<2G
Esequential/functional_1/batch_normalization_3/AssignMovingAvg_2/decay?
Nsequential/functional_1/batch_normalization_3/AssignMovingAvg_2/ReadVariableOpReadVariableOpEsequential_functional_1_batch_normalization_3_assignmovingavg_1225284R^sequential/functional_1/batch_normalization_3/AssignMovingAvg/AssignSubVariableOp*
_output_shapes
:*
dtype02P
Nsequential/functional_1/batch_normalization_3/AssignMovingAvg_2/ReadVariableOp?
Csequential/functional_1/batch_normalization_3/AssignMovingAvg_2/subSubVsequential/functional_1/batch_normalization_3/AssignMovingAvg_2/ReadVariableOp:value:0Hsequential/functional_1/batch_normalization_3/moments_1/Squeeze:output:0*
T0*X
_classN
LJloc:@sequential/functional_1/batch_normalization_3/AssignMovingAvg/1225284*
_output_shapes
:2E
Csequential/functional_1/batch_normalization_3/AssignMovingAvg_2/sub?
Csequential/functional_1/batch_normalization_3/AssignMovingAvg_2/mulMulGsequential/functional_1/batch_normalization_3/AssignMovingAvg_2/sub:z:0Nsequential/functional_1/batch_normalization_3/AssignMovingAvg_2/decay:output:0*
T0*X
_classN
LJloc:@sequential/functional_1/batch_normalization_3/AssignMovingAvg/1225284*
_output_shapes
:2E
Csequential/functional_1/batch_normalization_3/AssignMovingAvg_2/mul?
Ssequential/functional_1/batch_normalization_3/AssignMovingAvg_2/AssignSubVariableOpAssignSubVariableOpEsequential_functional_1_batch_normalization_3_assignmovingavg_1225284Gsequential/functional_1/batch_normalization_3/AssignMovingAvg_2/mul:z:0R^sequential/functional_1/batch_normalization_3/AssignMovingAvg/AssignSubVariableOpO^sequential/functional_1/batch_normalization_3/AssignMovingAvg_2/ReadVariableOp*X
_classN
LJloc:@sequential/functional_1/batch_normalization_3/AssignMovingAvg/1225284*
_output_shapes
 *
dtype02U
Ssequential/functional_1/batch_normalization_3/AssignMovingAvg_2/AssignSubVariableOp?
Esequential/functional_1/batch_normalization_3/AssignMovingAvg_3/decayConst*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_3/AssignMovingAvg_1/1225290*
_output_shapes
: *
dtype0*
valueB
 *
?#<2G
Esequential/functional_1/batch_normalization_3/AssignMovingAvg_3/decay?
Nsequential/functional_1/batch_normalization_3/AssignMovingAvg_3/ReadVariableOpReadVariableOpGsequential_functional_1_batch_normalization_3_assignmovingavg_1_1225290T^sequential/functional_1/batch_normalization_3/AssignMovingAvg_1/AssignSubVariableOp*
_output_shapes
:*
dtype02P
Nsequential/functional_1/batch_normalization_3/AssignMovingAvg_3/ReadVariableOp?
Csequential/functional_1/batch_normalization_3/AssignMovingAvg_3/subSubVsequential/functional_1/batch_normalization_3/AssignMovingAvg_3/ReadVariableOp:value:0Jsequential/functional_1/batch_normalization_3/moments_1/Squeeze_1:output:0*
T0*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_3/AssignMovingAvg_1/1225290*
_output_shapes
:2E
Csequential/functional_1/batch_normalization_3/AssignMovingAvg_3/sub?
Csequential/functional_1/batch_normalization_3/AssignMovingAvg_3/mulMulGsequential/functional_1/batch_normalization_3/AssignMovingAvg_3/sub:z:0Nsequential/functional_1/batch_normalization_3/AssignMovingAvg_3/decay:output:0*
T0*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_3/AssignMovingAvg_1/1225290*
_output_shapes
:2E
Csequential/functional_1/batch_normalization_3/AssignMovingAvg_3/mul?
Ssequential/functional_1/batch_normalization_3/AssignMovingAvg_3/AssignSubVariableOpAssignSubVariableOpGsequential_functional_1_batch_normalization_3_assignmovingavg_1_1225290Gsequential/functional_1/batch_normalization_3/AssignMovingAvg_3/mul:z:0T^sequential/functional_1/batch_normalization_3/AssignMovingAvg_1/AssignSubVariableOpO^sequential/functional_1/batch_normalization_3/AssignMovingAvg_3/ReadVariableOp*Z
_classP
NLloc:@sequential/functional_1/batch_normalization_3/AssignMovingAvg_1/1225290*
_output_shapes
 *
dtype02U
Ssequential/functional_1/batch_normalization_3/AssignMovingAvg_3/AssignSubVariableOp?
?sequential/functional_1/batch_normalization_3/batchnorm_1/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2A
?sequential/functional_1/batch_normalization_3/batchnorm_1/add/y?
=sequential/functional_1/batch_normalization_3/batchnorm_1/addAddV2Jsequential/functional_1/batch_normalization_3/moments_1/Squeeze_1:output:0Hsequential/functional_1/batch_normalization_3/batchnorm_1/add/y:output:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_3/batchnorm_1/add?
?sequential/functional_1/batch_normalization_3/batchnorm_1/RsqrtRsqrtAsequential/functional_1/batch_normalization_3/batchnorm_1/add:z:0*
T0*
_output_shapes
:2A
?sequential/functional_1/batch_normalization_3/batchnorm_1/Rsqrt?
Lsequential/functional_1/batch_normalization_3/batchnorm_1/mul/ReadVariableOpReadVariableOpSsequential_functional_1_batch_normalization_3_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02N
Lsequential/functional_1/batch_normalization_3/batchnorm_1/mul/ReadVariableOp?
=sequential/functional_1/batch_normalization_3/batchnorm_1/mulMulCsequential/functional_1/batch_normalization_3/batchnorm_1/Rsqrt:y:0Tsequential/functional_1/batch_normalization_3/batchnorm_1/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_3/batchnorm_1/mul?
?sequential/functional_1/batch_normalization_3/batchnorm_1/mul_1Mul4sequential/functional_1/conv3d_3/Elu_1:activations:0Asequential/functional_1/batch_normalization_3/batchnorm_1/mul:z:0*
T0*3
_output_shapes!
:?????????2A
?sequential/functional_1/batch_normalization_3/batchnorm_1/mul_1?
?sequential/functional_1/batch_normalization_3/batchnorm_1/mul_2MulHsequential/functional_1/batch_normalization_3/moments_1/Squeeze:output:0Asequential/functional_1/batch_normalization_3/batchnorm_1/mul:z:0*
T0*
_output_shapes
:2A
?sequential/functional_1/batch_normalization_3/batchnorm_1/mul_2?
Hsequential/functional_1/batch_normalization_3/batchnorm_1/ReadVariableOpReadVariableOpOsequential_functional_1_batch_normalization_3_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02J
Hsequential/functional_1/batch_normalization_3/batchnorm_1/ReadVariableOp?
=sequential/functional_1/batch_normalization_3/batchnorm_1/subSubPsequential/functional_1/batch_normalization_3/batchnorm_1/ReadVariableOp:value:0Csequential/functional_1/batch_normalization_3/batchnorm_1/mul_2:z:0*
T0*
_output_shapes
:2?
=sequential/functional_1/batch_normalization_3/batchnorm_1/sub?
?sequential/functional_1/batch_normalization_3/batchnorm_1/add_1AddV2Csequential/functional_1/batch_normalization_3/batchnorm_1/mul_1:z:0Asequential/functional_1/batch_normalization_3/batchnorm_1/sub:z:0*
T0*3
_output_shapes!
:?????????2A
?sequential/functional_1/batch_normalization_3/batchnorm_1/add_1?
5sequential/functional_1/average_pooling3d/AvgPool3D_1	AvgPool3DCsequential/functional_1/batch_normalization_3/batchnorm_1/add_1:z:0*
T0*3
_output_shapes!
:?????????*
ksize	
*
paddingVALID*
strides	
27
5sequential/functional_1/average_pooling3d/AvgPool3D_1?
'sequential/functional_1/flatten/Const_1Const*
_output_shapes
:*
dtype0*
valueB"????   2)
'sequential/functional_1/flatten/Const_1?
)sequential/functional_1/flatten/Reshape_1Reshape>sequential/functional_1/average_pooling3d/AvgPool3D_1:output:00sequential/functional_1/flatten/Const_1:output:0*
T0*(
_output_shapes
:??????????
2+
)sequential/functional_1/flatten/Reshape_1?
/sequential/functional_1/dropout/dropout_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *UU??21
/sequential/functional_1/dropout/dropout_1/Const?
-sequential/functional_1/dropout/dropout_1/MulMul2sequential/functional_1/flatten/Reshape_1:output:08sequential/functional_1/dropout/dropout_1/Const:output:0*
T0*(
_output_shapes
:??????????
2/
-sequential/functional_1/dropout/dropout_1/Mul?
/sequential/functional_1/dropout/dropout_1/ShapeShape2sequential/functional_1/flatten/Reshape_1:output:0*
T0*
_output_shapes
:21
/sequential/functional_1/dropout/dropout_1/Shape?
Fsequential/functional_1/dropout/dropout_1/random_uniform/RandomUniformRandomUniform8sequential/functional_1/dropout/dropout_1/Shape:output:0*
T0*(
_output_shapes
:??????????
*
dtype02H
Fsequential/functional_1/dropout/dropout_1/random_uniform/RandomUniform?
8sequential/functional_1/dropout/dropout_1/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2:
8sequential/functional_1/dropout/dropout_1/GreaterEqual/y?
6sequential/functional_1/dropout/dropout_1/GreaterEqualGreaterEqualOsequential/functional_1/dropout/dropout_1/random_uniform/RandomUniform:output:0Asequential/functional_1/dropout/dropout_1/GreaterEqual/y:output:0*
T0*(
_output_shapes
:??????????
28
6sequential/functional_1/dropout/dropout_1/GreaterEqual?
.sequential/functional_1/dropout/dropout_1/CastCast:sequential/functional_1/dropout/dropout_1/GreaterEqual:z:0*

DstT0*

SrcT0
*(
_output_shapes
:??????????
20
.sequential/functional_1/dropout/dropout_1/Cast?
/sequential/functional_1/dropout/dropout_1/Mul_1Mul1sequential/functional_1/dropout/dropout_1/Mul:z:02sequential/functional_1/dropout/dropout_1/Cast:y:0*
T0*(
_output_shapes
:??????????
21
/sequential/functional_1/dropout/dropout_1/Mul_1?
6sequential/functional_1/layer1/MatMul_1/ReadVariableOpReadVariableOp=sequential_functional_1_layer1_matmul_readvariableop_resource* 
_output_shapes
:
?
?*
dtype028
6sequential/functional_1/layer1/MatMul_1/ReadVariableOp?
'sequential/functional_1/layer1/MatMul_1MatMul3sequential/functional_1/dropout/dropout_1/Mul_1:z:0>sequential/functional_1/layer1/MatMul_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2)
'sequential/functional_1/layer1/MatMul_1?
7sequential/functional_1/layer1/BiasAdd_1/ReadVariableOpReadVariableOp>sequential_functional_1_layer1_biasadd_readvariableop_resource*
_output_shapes	
:?*
dtype029
7sequential/functional_1/layer1/BiasAdd_1/ReadVariableOp?
(sequential/functional_1/layer1/BiasAdd_1BiasAdd1sequential/functional_1/layer1/MatMul_1:product:0?sequential/functional_1/layer1/BiasAdd_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2*
(sequential/functional_1/layer1/BiasAdd_1?
$sequential/functional_1/layer1/Elu_1Elu1sequential/functional_1/layer1/BiasAdd_1:output:0*
T0*(
_output_shapes
:??????????2&
$sequential/functional_1/layer1/Elu_1?
(sequential/dense/MatMul_1/ReadVariableOpReadVariableOp/sequential_dense_matmul_readvariableop_resource*
_output_shapes
:	?d*
dtype02*
(sequential/dense/MatMul_1/ReadVariableOp?
sequential/dense/MatMul_1MatMul2sequential/functional_1/layer1/Elu_1:activations:00sequential/dense/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????d2
sequential/dense/MatMul_1?
)sequential/dense/BiasAdd_1/ReadVariableOpReadVariableOp0sequential_dense_biasadd_readvariableop_resource*
_output_shapes
:d*
dtype02+
)sequential/dense/BiasAdd_1/ReadVariableOp?
sequential/dense/BiasAdd_1BiasAdd#sequential/dense/MatMul_1:product:01sequential/dense/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????d2
sequential/dense/BiasAdd_1?
sequential/dense/Elu_1Elu#sequential/dense/BiasAdd_1:output:0*
T0*'
_output_shapes
:?????????d2
sequential/dense/Elu_1?
*sequential/dense_1/MatMul_1/ReadVariableOpReadVariableOp1sequential_dense_1_matmul_readvariableop_resource*
_output_shapes

:d
*
dtype02,
*sequential/dense_1/MatMul_1/ReadVariableOp?
sequential/dense_1/MatMul_1MatMul$sequential/dense/Elu_1:activations:02sequential/dense_1/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????
2
sequential/dense_1/MatMul_1?
+sequential/dense_1/BiasAdd_1/ReadVariableOpReadVariableOp2sequential_dense_1_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02-
+sequential/dense_1/BiasAdd_1/ReadVariableOp?
sequential/dense_1/BiasAdd_1BiasAdd%sequential/dense_1/MatMul_1:product:03sequential/dense_1/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????
2
sequential/dense_1/BiasAdd_1?
sequential/dense_1/Elu_1Elu%sequential/dense_1/BiasAdd_1:output:0*
T0*'
_output_shapes
:?????????
2
sequential/dense_1/Elu_1?

lambda/subSub$sequential/dense_1/Elu:activations:0&sequential/dense_1/Elu_1:activations:0*
T0*'
_output_shapes
:?????????
2

lambda/suba

lambda/AbsAbslambda/sub:z:0*
T0*'
_output_shapes
:?????????
2

lambda/Abst
concatenate/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
concatenate/concat/axis?
concatenate/concatConcatV2lambda/Abs:y:0inputs_2 concatenate/concat/axis:output:0*
N*
T0*'
_output_shapes
:?????????2
concatenate/concat?
dense_2/MatMul/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes

:*
dtype02
dense_2/MatMul/ReadVariableOp?
dense_2/MatMulMatMulconcatenate/concat:output:0%dense_2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
dense_2/MatMul?
dense_2/BiasAdd/ReadVariableOpReadVariableOp'dense_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02 
dense_2/BiasAdd/ReadVariableOp?
dense_2/BiasAddBiasAdddense_2/MatMul:product:0&dense_2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
dense_2/BiasAdd?
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOp/sequential_dense_matmul_readvariableop_resource*
_output_shapes
:	?d*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOp?
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?d2!
dense/kernel/Regularizer/Square?
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const?
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/Sum?
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2 
dense/kernel/Regularizer/mul/x?
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/mul?
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp1sequential_dense_1_matmul_readvariableop_resource*
_output_shapes

:d
*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp?
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d
2#
!dense_1/kernel/Regularizer/Square?
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const?
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/Sum?
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2"
 dense_1/kernel/Regularizer/mul/x?
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/mul?
IdentityIdentitydense_2/BiasAdd:output:0P^sequential/functional_1/batch_normalization/AssignMovingAvg/AssignSubVariableOpR^sequential/functional_1/batch_normalization/AssignMovingAvg_1/AssignSubVariableOpR^sequential/functional_1/batch_normalization/AssignMovingAvg_2/AssignSubVariableOpR^sequential/functional_1/batch_normalization/AssignMovingAvg_3/AssignSubVariableOpR^sequential/functional_1/batch_normalization_1/AssignMovingAvg/AssignSubVariableOpT^sequential/functional_1/batch_normalization_1/AssignMovingAvg_1/AssignSubVariableOpT^sequential/functional_1/batch_normalization_1/AssignMovingAvg_2/AssignSubVariableOpT^sequential/functional_1/batch_normalization_1/AssignMovingAvg_3/AssignSubVariableOpR^sequential/functional_1/batch_normalization_2/AssignMovingAvg/AssignSubVariableOpT^sequential/functional_1/batch_normalization_2/AssignMovingAvg_1/AssignSubVariableOpT^sequential/functional_1/batch_normalization_2/AssignMovingAvg_2/AssignSubVariableOpT^sequential/functional_1/batch_normalization_2/AssignMovingAvg_3/AssignSubVariableOpR^sequential/functional_1/batch_normalization_3/AssignMovingAvg/AssignSubVariableOpT^sequential/functional_1/batch_normalization_3/AssignMovingAvg_1/AssignSubVariableOpT^sequential/functional_1/batch_normalization_3/AssignMovingAvg_2/AssignSubVariableOpT^sequential/functional_1/batch_normalization_3/AssignMovingAvg_3/AssignSubVariableOp*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????:??????????:?????????
::::::::::::::::::::::::::::::::2?
Osequential/functional_1/batch_normalization/AssignMovingAvg/AssignSubVariableOpOsequential/functional_1/batch_normalization/AssignMovingAvg/AssignSubVariableOp2?
Qsequential/functional_1/batch_normalization/AssignMovingAvg_1/AssignSubVariableOpQsequential/functional_1/batch_normalization/AssignMovingAvg_1/AssignSubVariableOp2?
Qsequential/functional_1/batch_normalization/AssignMovingAvg_2/AssignSubVariableOpQsequential/functional_1/batch_normalization/AssignMovingAvg_2/AssignSubVariableOp2?
Qsequential/functional_1/batch_normalization/AssignMovingAvg_3/AssignSubVariableOpQsequential/functional_1/batch_normalization/AssignMovingAvg_3/AssignSubVariableOp2?
Qsequential/functional_1/batch_normalization_1/AssignMovingAvg/AssignSubVariableOpQsequential/functional_1/batch_normalization_1/AssignMovingAvg/AssignSubVariableOp2?
Ssequential/functional_1/batch_normalization_1/AssignMovingAvg_1/AssignSubVariableOpSsequential/functional_1/batch_normalization_1/AssignMovingAvg_1/AssignSubVariableOp2?
Ssequential/functional_1/batch_normalization_1/AssignMovingAvg_2/AssignSubVariableOpSsequential/functional_1/batch_normalization_1/AssignMovingAvg_2/AssignSubVariableOp2?
Ssequential/functional_1/batch_normalization_1/AssignMovingAvg_3/AssignSubVariableOpSsequential/functional_1/batch_normalization_1/AssignMovingAvg_3/AssignSubVariableOp2?
Qsequential/functional_1/batch_normalization_2/AssignMovingAvg/AssignSubVariableOpQsequential/functional_1/batch_normalization_2/AssignMovingAvg/AssignSubVariableOp2?
Ssequential/functional_1/batch_normalization_2/AssignMovingAvg_1/AssignSubVariableOpSsequential/functional_1/batch_normalization_2/AssignMovingAvg_1/AssignSubVariableOp2?
Ssequential/functional_1/batch_normalization_2/AssignMovingAvg_2/AssignSubVariableOpSsequential/functional_1/batch_normalization_2/AssignMovingAvg_2/AssignSubVariableOp2?
Ssequential/functional_1/batch_normalization_2/AssignMovingAvg_3/AssignSubVariableOpSsequential/functional_1/batch_normalization_2/AssignMovingAvg_3/AssignSubVariableOp2?
Qsequential/functional_1/batch_normalization_3/AssignMovingAvg/AssignSubVariableOpQsequential/functional_1/batch_normalization_3/AssignMovingAvg/AssignSubVariableOp2?
Ssequential/functional_1/batch_normalization_3/AssignMovingAvg_1/AssignSubVariableOpSsequential/functional_1/batch_normalization_3/AssignMovingAvg_1/AssignSubVariableOp2?
Ssequential/functional_1/batch_normalization_3/AssignMovingAvg_2/AssignSubVariableOpSsequential/functional_1/batch_normalization_3/AssignMovingAvg_2/AssignSubVariableOp2?
Ssequential/functional_1/batch_normalization_3/AssignMovingAvg_3/AssignSubVariableOpSsequential/functional_1/batch_normalization_3/AssignMovingAvg_3/AssignSubVariableOp:^ Z
4
_output_shapes"
 :??????????
"
_user_specified_name
inputs/0:^Z
4
_output_shapes"
 :??????????
"
_user_specified_name
inputs/1:QM
'
_output_shapes
:?????????

"
_user_specified_name
inputs/2
?
`
D__inference_flatten_layer_call_and_return_conditional_losses_1227651

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"????   2
Consth
ReshapeReshapeinputsConst:output:0*
T0*(
_output_shapes
:??????????
2	
Reshapee
IdentityIdentityReshape:output:0*
T0*(
_output_shapes
:??????????
2

Identity"
identityIdentity:output:0*2
_input_shapes!
:?????????:[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?*
?
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_1223065

inputs
assignmovingavg_1223040
assignmovingavg_1_1223046)
%batchnorm_mul_readvariableop_resource%
!batchnorm_readvariableop_resource
identity??#AssignMovingAvg/AssignSubVariableOp?%AssignMovingAvg_1/AssignSubVariableOp?
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2 
moments/mean/reduction_indices?
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/mean?
moments/StopGradientStopGradientmoments/mean:output:0*
T0**
_output_shapes
:2
moments/StopGradient?
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*3
_output_shapes!
:?????????2
moments/SquaredDifference?
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2$
"moments/variance/reduction_indices?
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/variance?
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze?
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze_1?
AssignMovingAvg/decayConst**
_class 
loc:@AssignMovingAvg/1223040*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg/decay?
AssignMovingAvg/ReadVariableOpReadVariableOpassignmovingavg_1223040*
_output_shapes
:*
dtype02 
AssignMovingAvg/ReadVariableOp?
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0**
_class 
loc:@AssignMovingAvg/1223040*
_output_shapes
:2
AssignMovingAvg/sub?
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0**
_class 
loc:@AssignMovingAvg/1223040*
_output_shapes
:2
AssignMovingAvg/mul?
#AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1223040AssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp**
_class 
loc:@AssignMovingAvg/1223040*
_output_shapes
 *
dtype02%
#AssignMovingAvg/AssignSubVariableOp?
AssignMovingAvg_1/decayConst*,
_class"
 loc:@AssignMovingAvg_1/1223046*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg_1/decay?
 AssignMovingAvg_1/ReadVariableOpReadVariableOpassignmovingavg_1_1223046*
_output_shapes
:*
dtype02"
 AssignMovingAvg_1/ReadVariableOp?
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1223046*
_output_shapes
:2
AssignMovingAvg_1/sub?
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1223046*
_output_shapes
:2
AssignMovingAvg_1/mul?
%AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1_1223046AssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*,
_class"
 loc:@AssignMovingAvg_1/1223046*
_output_shapes
 *
dtype02'
%AssignMovingAvg_1/AssignSubVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/mul_1{
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp?
batchnorm/subSub batchnorm/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0$^AssignMovingAvg/AssignSubVariableOp&^AssignMovingAvg_1/AssignSubVariableOp*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????::::2J
#AssignMovingAvg/AssignSubVariableOp#AssignMovingAvg/AssignSubVariableOp2N
%AssignMovingAvg_1/AssignSubVariableOp%AssignMovingAvg_1/AssignSubVariableOp:[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?	
?
E__inference_conv3d_2_layer_call_and_return_conditional_losses_1227288

inputs"
conv3d_readvariableop_resource#
biasadd_readvariableop_resource
identity??
Conv3D/ReadVariableOpReadVariableOpconv3d_readvariableop_resource**
_output_shapes
:*
dtype02
Conv3D/ReadVariableOp?
Conv3DConv3DinputsConv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
Conv3D?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddConv3D:output:0BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2	
BiasAdda
EluEluBiasAdd:output:0*
T0*3
_output_shapes!
:?????????2
Eluq
IdentityIdentityElu:activations:0*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*:
_input_shapes)
':?????????:::[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?*
?
P__inference_batch_normalization_layer_call_and_return_conditional_losses_1222829

inputs
assignmovingavg_1222804
assignmovingavg_1_1222810)
%batchnorm_mul_readvariableop_resource%
!batchnorm_readvariableop_resource
identity??#AssignMovingAvg/AssignSubVariableOp?%AssignMovingAvg_1/AssignSubVariableOp?
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2 
moments/mean/reduction_indices?
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/mean?
moments/StopGradientStopGradientmoments/mean:output:0*
T0**
_output_shapes
:2
moments/StopGradient?
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*3
_output_shapes!
:?????????2
moments/SquaredDifference?
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2$
"moments/variance/reduction_indices?
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/variance?
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze?
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze_1?
AssignMovingAvg/decayConst**
_class 
loc:@AssignMovingAvg/1222804*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg/decay?
AssignMovingAvg/ReadVariableOpReadVariableOpassignmovingavg_1222804*
_output_shapes
:*
dtype02 
AssignMovingAvg/ReadVariableOp?
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0**
_class 
loc:@AssignMovingAvg/1222804*
_output_shapes
:2
AssignMovingAvg/sub?
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0**
_class 
loc:@AssignMovingAvg/1222804*
_output_shapes
:2
AssignMovingAvg/mul?
#AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1222804AssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp**
_class 
loc:@AssignMovingAvg/1222804*
_output_shapes
 *
dtype02%
#AssignMovingAvg/AssignSubVariableOp?
AssignMovingAvg_1/decayConst*,
_class"
 loc:@AssignMovingAvg_1/1222810*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg_1/decay?
 AssignMovingAvg_1/ReadVariableOpReadVariableOpassignmovingavg_1_1222810*
_output_shapes
:*
dtype02"
 AssignMovingAvg_1/ReadVariableOp?
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1222810*
_output_shapes
:2
AssignMovingAvg_1/sub?
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1222810*
_output_shapes
:2
AssignMovingAvg_1/mul?
%AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1_1222810AssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*,
_class"
 loc:@AssignMovingAvg_1/1222810*
_output_shapes
 *
dtype02'
%AssignMovingAvg_1/AssignSubVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/mul_1{
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp?
batchnorm/subSub batchnorm/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0$^AssignMovingAvg/AssignSubVariableOp&^AssignMovingAvg_1/AssignSubVariableOp*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????::::2J
#AssignMovingAvg/AssignSubVariableOp#AssignMovingAvg/AssignSubVariableOp2N
%AssignMovingAvg_1/AssignSubVariableOp%AssignMovingAvg_1/AssignSubVariableOp:[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?
Y
-__inference_concatenate_layer_call_fn_1226408
inputs_0
inputs_1
identity?
PartitionedCallPartitionedCallinputs_0inputs_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *Q
fLRJ
H__inference_concatenate_layer_call_and_return_conditional_losses_12245132
PartitionedCalll
IdentityIdentityPartitionedCall:output:0*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*9
_input_shapes(
&:?????????
:?????????
:Q M
'
_output_shapes
:?????????

"
_user_specified_name
inputs/0:QM
'
_output_shapes
:?????????

"
_user_specified_name
inputs/1
??
?
"__inference__wrapped_model_1222192
input_1
input_2
input_3N
Jfunctional_3_sequential_functional_1_conv3d_conv3d_readvariableop_resourceO
Kfunctional_3_sequential_functional_1_conv3d_biasadd_readvariableop_resource^
Zfunctional_3_sequential_functional_1_batch_normalization_batchnorm_readvariableop_resourceb
^functional_3_sequential_functional_1_batch_normalization_batchnorm_mul_readvariableop_resource`
\functional_3_sequential_functional_1_batch_normalization_batchnorm_readvariableop_1_resource`
\functional_3_sequential_functional_1_batch_normalization_batchnorm_readvariableop_2_resourceP
Lfunctional_3_sequential_functional_1_conv3d_1_conv3d_readvariableop_resourceQ
Mfunctional_3_sequential_functional_1_conv3d_1_biasadd_readvariableop_resource`
\functional_3_sequential_functional_1_batch_normalization_1_batchnorm_readvariableop_resourced
`functional_3_sequential_functional_1_batch_normalization_1_batchnorm_mul_readvariableop_resourceb
^functional_3_sequential_functional_1_batch_normalization_1_batchnorm_readvariableop_1_resourceb
^functional_3_sequential_functional_1_batch_normalization_1_batchnorm_readvariableop_2_resourceP
Lfunctional_3_sequential_functional_1_conv3d_2_conv3d_readvariableop_resourceQ
Mfunctional_3_sequential_functional_1_conv3d_2_biasadd_readvariableop_resource`
\functional_3_sequential_functional_1_batch_normalization_2_batchnorm_readvariableop_resourced
`functional_3_sequential_functional_1_batch_normalization_2_batchnorm_mul_readvariableop_resourceb
^functional_3_sequential_functional_1_batch_normalization_2_batchnorm_readvariableop_1_resourceb
^functional_3_sequential_functional_1_batch_normalization_2_batchnorm_readvariableop_2_resourceP
Lfunctional_3_sequential_functional_1_conv3d_3_conv3d_readvariableop_resourceQ
Mfunctional_3_sequential_functional_1_conv3d_3_biasadd_readvariableop_resource`
\functional_3_sequential_functional_1_batch_normalization_3_batchnorm_readvariableop_resourced
`functional_3_sequential_functional_1_batch_normalization_3_batchnorm_mul_readvariableop_resourceb
^functional_3_sequential_functional_1_batch_normalization_3_batchnorm_readvariableop_1_resourceb
^functional_3_sequential_functional_1_batch_normalization_3_batchnorm_readvariableop_2_resourceN
Jfunctional_3_sequential_functional_1_layer1_matmul_readvariableop_resourceO
Kfunctional_3_sequential_functional_1_layer1_biasadd_readvariableop_resource@
<functional_3_sequential_dense_matmul_readvariableop_resourceA
=functional_3_sequential_dense_biasadd_readvariableop_resourceB
>functional_3_sequential_dense_1_matmul_readvariableop_resourceC
?functional_3_sequential_dense_1_biasadd_readvariableop_resource7
3functional_3_dense_2_matmul_readvariableop_resource8
4functional_3_dense_2_biasadd_readvariableop_resource
identity??
Afunctional_3/sequential/functional_1/conv3d/Conv3D/ReadVariableOpReadVariableOpJfunctional_3_sequential_functional_1_conv3d_conv3d_readvariableop_resource*+
_output_shapes
:?*
dtype02C
Afunctional_3/sequential/functional_1/conv3d/Conv3D/ReadVariableOp?
2functional_3/sequential/functional_1/conv3d/Conv3DConv3Dinput_1Ifunctional_3/sequential/functional_1/conv3d/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
24
2functional_3/sequential/functional_1/conv3d/Conv3D?
Bfunctional_3/sequential/functional_1/conv3d/BiasAdd/ReadVariableOpReadVariableOpKfunctional_3_sequential_functional_1_conv3d_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02D
Bfunctional_3/sequential/functional_1/conv3d/BiasAdd/ReadVariableOp?
3functional_3/sequential/functional_1/conv3d/BiasAddBiasAdd;functional_3/sequential/functional_1/conv3d/Conv3D:output:0Jfunctional_3/sequential/functional_1/conv3d/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????25
3functional_3/sequential/functional_1/conv3d/BiasAdd?
Qfunctional_3/sequential/functional_1/batch_normalization/batchnorm/ReadVariableOpReadVariableOpZfunctional_3_sequential_functional_1_batch_normalization_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02S
Qfunctional_3/sequential/functional_1/batch_normalization/batchnorm/ReadVariableOp?
Hfunctional_3/sequential/functional_1/batch_normalization/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2J
Hfunctional_3/sequential/functional_1/batch_normalization/batchnorm/add/y?
Ffunctional_3/sequential/functional_1/batch_normalization/batchnorm/addAddV2Yfunctional_3/sequential/functional_1/batch_normalization/batchnorm/ReadVariableOp:value:0Qfunctional_3/sequential/functional_1/batch_normalization/batchnorm/add/y:output:0*
T0*
_output_shapes
:2H
Ffunctional_3/sequential/functional_1/batch_normalization/batchnorm/add?
Hfunctional_3/sequential/functional_1/batch_normalization/batchnorm/RsqrtRsqrtJfunctional_3/sequential/functional_1/batch_normalization/batchnorm/add:z:0*
T0*
_output_shapes
:2J
Hfunctional_3/sequential/functional_1/batch_normalization/batchnorm/Rsqrt?
Ufunctional_3/sequential/functional_1/batch_normalization/batchnorm/mul/ReadVariableOpReadVariableOp^functional_3_sequential_functional_1_batch_normalization_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02W
Ufunctional_3/sequential/functional_1/batch_normalization/batchnorm/mul/ReadVariableOp?
Ffunctional_3/sequential/functional_1/batch_normalization/batchnorm/mulMulLfunctional_3/sequential/functional_1/batch_normalization/batchnorm/Rsqrt:y:0]functional_3/sequential/functional_1/batch_normalization/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2H
Ffunctional_3/sequential/functional_1/batch_normalization/batchnorm/mul?
Hfunctional_3/sequential/functional_1/batch_normalization/batchnorm/mul_1Mul<functional_3/sequential/functional_1/conv3d/BiasAdd:output:0Jfunctional_3/sequential/functional_1/batch_normalization/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2J
Hfunctional_3/sequential/functional_1/batch_normalization/batchnorm/mul_1?
Sfunctional_3/sequential/functional_1/batch_normalization/batchnorm/ReadVariableOp_1ReadVariableOp\functional_3_sequential_functional_1_batch_normalization_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02U
Sfunctional_3/sequential/functional_1/batch_normalization/batchnorm/ReadVariableOp_1?
Hfunctional_3/sequential/functional_1/batch_normalization/batchnorm/mul_2Mul[functional_3/sequential/functional_1/batch_normalization/batchnorm/ReadVariableOp_1:value:0Jfunctional_3/sequential/functional_1/batch_normalization/batchnorm/mul:z:0*
T0*
_output_shapes
:2J
Hfunctional_3/sequential/functional_1/batch_normalization/batchnorm/mul_2?
Sfunctional_3/sequential/functional_1/batch_normalization/batchnorm/ReadVariableOp_2ReadVariableOp\functional_3_sequential_functional_1_batch_normalization_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02U
Sfunctional_3/sequential/functional_1/batch_normalization/batchnorm/ReadVariableOp_2?
Ffunctional_3/sequential/functional_1/batch_normalization/batchnorm/subSub[functional_3/sequential/functional_1/batch_normalization/batchnorm/ReadVariableOp_2:value:0Lfunctional_3/sequential/functional_1/batch_normalization/batchnorm/mul_2:z:0*
T0*
_output_shapes
:2H
Ffunctional_3/sequential/functional_1/batch_normalization/batchnorm/sub?
Hfunctional_3/sequential/functional_1/batch_normalization/batchnorm/add_1AddV2Lfunctional_3/sequential/functional_1/batch_normalization/batchnorm/mul_1:z:0Jfunctional_3/sequential/functional_1/batch_normalization/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2J
Hfunctional_3/sequential/functional_1/batch_normalization/batchnorm/add_1?
Cfunctional_3/sequential/functional_1/conv3d_1/Conv3D/ReadVariableOpReadVariableOpLfunctional_3_sequential_functional_1_conv3d_1_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02E
Cfunctional_3/sequential/functional_1/conv3d_1/Conv3D/ReadVariableOp?
4functional_3/sequential/functional_1/conv3d_1/Conv3DConv3DLfunctional_3/sequential/functional_1/batch_normalization/batchnorm/add_1:z:0Kfunctional_3/sequential/functional_1/conv3d_1/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
26
4functional_3/sequential/functional_1/conv3d_1/Conv3D?
Dfunctional_3/sequential/functional_1/conv3d_1/BiasAdd/ReadVariableOpReadVariableOpMfunctional_3_sequential_functional_1_conv3d_1_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02F
Dfunctional_3/sequential/functional_1/conv3d_1/BiasAdd/ReadVariableOp?
5functional_3/sequential/functional_1/conv3d_1/BiasAddBiasAdd=functional_3/sequential/functional_1/conv3d_1/Conv3D:output:0Lfunctional_3/sequential/functional_1/conv3d_1/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????27
5functional_3/sequential/functional_1/conv3d_1/BiasAdd?
1functional_3/sequential/functional_1/conv3d_1/EluElu>functional_3/sequential/functional_1/conv3d_1/BiasAdd:output:0*
T0*3
_output_shapes!
:?????????23
1functional_3/sequential/functional_1/conv3d_1/Elu?
Sfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/ReadVariableOpReadVariableOp\functional_3_sequential_functional_1_batch_normalization_1_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02U
Sfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/ReadVariableOp?
Jfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2L
Jfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/add/y?
Hfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/addAddV2[functional_3/sequential/functional_1/batch_normalization_1/batchnorm/ReadVariableOp:value:0Sfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/add/y:output:0*
T0*
_output_shapes
:2J
Hfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/add?
Jfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/RsqrtRsqrtLfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/add:z:0*
T0*
_output_shapes
:2L
Jfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/Rsqrt?
Wfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/mul/ReadVariableOpReadVariableOp`functional_3_sequential_functional_1_batch_normalization_1_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02Y
Wfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/mul/ReadVariableOp?
Hfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/mulMulNfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/Rsqrt:y:0_functional_3/sequential/functional_1/batch_normalization_1/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2J
Hfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/mul?
Jfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/mul_1Mul?functional_3/sequential/functional_1/conv3d_1/Elu:activations:0Lfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2L
Jfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/mul_1?
Ufunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/ReadVariableOp_1ReadVariableOp^functional_3_sequential_functional_1_batch_normalization_1_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02W
Ufunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/ReadVariableOp_1?
Jfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/mul_2Mul]functional_3/sequential/functional_1/batch_normalization_1/batchnorm/ReadVariableOp_1:value:0Lfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/mul:z:0*
T0*
_output_shapes
:2L
Jfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/mul_2?
Ufunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/ReadVariableOp_2ReadVariableOp^functional_3_sequential_functional_1_batch_normalization_1_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02W
Ufunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/ReadVariableOp_2?
Hfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/subSub]functional_3/sequential/functional_1/batch_normalization_1/batchnorm/ReadVariableOp_2:value:0Nfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/mul_2:z:0*
T0*
_output_shapes
:2J
Hfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/sub?
Jfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/add_1AddV2Nfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/mul_1:z:0Lfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2L
Jfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/add_1?
Cfunctional_3/sequential/functional_1/conv3d_2/Conv3D/ReadVariableOpReadVariableOpLfunctional_3_sequential_functional_1_conv3d_2_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02E
Cfunctional_3/sequential/functional_1/conv3d_2/Conv3D/ReadVariableOp?
4functional_3/sequential/functional_1/conv3d_2/Conv3DConv3DNfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm/add_1:z:0Kfunctional_3/sequential/functional_1/conv3d_2/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
26
4functional_3/sequential/functional_1/conv3d_2/Conv3D?
Dfunctional_3/sequential/functional_1/conv3d_2/BiasAdd/ReadVariableOpReadVariableOpMfunctional_3_sequential_functional_1_conv3d_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02F
Dfunctional_3/sequential/functional_1/conv3d_2/BiasAdd/ReadVariableOp?
5functional_3/sequential/functional_1/conv3d_2/BiasAddBiasAdd=functional_3/sequential/functional_1/conv3d_2/Conv3D:output:0Lfunctional_3/sequential/functional_1/conv3d_2/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????27
5functional_3/sequential/functional_1/conv3d_2/BiasAdd?
1functional_3/sequential/functional_1/conv3d_2/EluElu>functional_3/sequential/functional_1/conv3d_2/BiasAdd:output:0*
T0*3
_output_shapes!
:?????????23
1functional_3/sequential/functional_1/conv3d_2/Elu?
Sfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/ReadVariableOpReadVariableOp\functional_3_sequential_functional_1_batch_normalization_2_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02U
Sfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/ReadVariableOp?
Jfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2L
Jfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/add/y?
Hfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/addAddV2[functional_3/sequential/functional_1/batch_normalization_2/batchnorm/ReadVariableOp:value:0Sfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/add/y:output:0*
T0*
_output_shapes
:2J
Hfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/add?
Jfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/RsqrtRsqrtLfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/add:z:0*
T0*
_output_shapes
:2L
Jfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/Rsqrt?
Wfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/mul/ReadVariableOpReadVariableOp`functional_3_sequential_functional_1_batch_normalization_2_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02Y
Wfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/mul/ReadVariableOp?
Hfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/mulMulNfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/Rsqrt:y:0_functional_3/sequential/functional_1/batch_normalization_2/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2J
Hfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/mul?
Jfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/mul_1Mul?functional_3/sequential/functional_1/conv3d_2/Elu:activations:0Lfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2L
Jfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/mul_1?
Ufunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/ReadVariableOp_1ReadVariableOp^functional_3_sequential_functional_1_batch_normalization_2_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02W
Ufunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/ReadVariableOp_1?
Jfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/mul_2Mul]functional_3/sequential/functional_1/batch_normalization_2/batchnorm/ReadVariableOp_1:value:0Lfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/mul:z:0*
T0*
_output_shapes
:2L
Jfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/mul_2?
Ufunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/ReadVariableOp_2ReadVariableOp^functional_3_sequential_functional_1_batch_normalization_2_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02W
Ufunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/ReadVariableOp_2?
Hfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/subSub]functional_3/sequential/functional_1/batch_normalization_2/batchnorm/ReadVariableOp_2:value:0Nfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/mul_2:z:0*
T0*
_output_shapes
:2J
Hfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/sub?
Jfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/add_1AddV2Nfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/mul_1:z:0Lfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2L
Jfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/add_1?
Cfunctional_3/sequential/functional_1/conv3d_3/Conv3D/ReadVariableOpReadVariableOpLfunctional_3_sequential_functional_1_conv3d_3_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02E
Cfunctional_3/sequential/functional_1/conv3d_3/Conv3D/ReadVariableOp?
4functional_3/sequential/functional_1/conv3d_3/Conv3DConv3DNfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm/add_1:z:0Kfunctional_3/sequential/functional_1/conv3d_3/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
26
4functional_3/sequential/functional_1/conv3d_3/Conv3D?
Dfunctional_3/sequential/functional_1/conv3d_3/BiasAdd/ReadVariableOpReadVariableOpMfunctional_3_sequential_functional_1_conv3d_3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02F
Dfunctional_3/sequential/functional_1/conv3d_3/BiasAdd/ReadVariableOp?
5functional_3/sequential/functional_1/conv3d_3/BiasAddBiasAdd=functional_3/sequential/functional_1/conv3d_3/Conv3D:output:0Lfunctional_3/sequential/functional_1/conv3d_3/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????27
5functional_3/sequential/functional_1/conv3d_3/BiasAdd?
1functional_3/sequential/functional_1/conv3d_3/EluElu>functional_3/sequential/functional_1/conv3d_3/BiasAdd:output:0*
T0*3
_output_shapes!
:?????????23
1functional_3/sequential/functional_1/conv3d_3/Elu?
Sfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/ReadVariableOpReadVariableOp\functional_3_sequential_functional_1_batch_normalization_3_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02U
Sfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/ReadVariableOp?
Jfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2L
Jfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/add/y?
Hfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/addAddV2[functional_3/sequential/functional_1/batch_normalization_3/batchnorm/ReadVariableOp:value:0Sfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/add/y:output:0*
T0*
_output_shapes
:2J
Hfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/add?
Jfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/RsqrtRsqrtLfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/add:z:0*
T0*
_output_shapes
:2L
Jfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/Rsqrt?
Wfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/mul/ReadVariableOpReadVariableOp`functional_3_sequential_functional_1_batch_normalization_3_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02Y
Wfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/mul/ReadVariableOp?
Hfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/mulMulNfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/Rsqrt:y:0_functional_3/sequential/functional_1/batch_normalization_3/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2J
Hfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/mul?
Jfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/mul_1Mul?functional_3/sequential/functional_1/conv3d_3/Elu:activations:0Lfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2L
Jfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/mul_1?
Ufunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/ReadVariableOp_1ReadVariableOp^functional_3_sequential_functional_1_batch_normalization_3_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02W
Ufunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/ReadVariableOp_1?
Jfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/mul_2Mul]functional_3/sequential/functional_1/batch_normalization_3/batchnorm/ReadVariableOp_1:value:0Lfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/mul:z:0*
T0*
_output_shapes
:2L
Jfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/mul_2?
Ufunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/ReadVariableOp_2ReadVariableOp^functional_3_sequential_functional_1_batch_normalization_3_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02W
Ufunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/ReadVariableOp_2?
Hfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/subSub]functional_3/sequential/functional_1/batch_normalization_3/batchnorm/ReadVariableOp_2:value:0Nfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/mul_2:z:0*
T0*
_output_shapes
:2J
Hfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/sub?
Jfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/add_1AddV2Nfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/mul_1:z:0Lfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2L
Jfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/add_1?
@functional_3/sequential/functional_1/average_pooling3d/AvgPool3D	AvgPool3DNfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm/add_1:z:0*
T0*3
_output_shapes!
:?????????*
ksize	
*
paddingVALID*
strides	
2B
@functional_3/sequential/functional_1/average_pooling3d/AvgPool3D?
2functional_3/sequential/functional_1/flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"????   24
2functional_3/sequential/functional_1/flatten/Const?
4functional_3/sequential/functional_1/flatten/ReshapeReshapeIfunctional_3/sequential/functional_1/average_pooling3d/AvgPool3D:output:0;functional_3/sequential/functional_1/flatten/Const:output:0*
T0*(
_output_shapes
:??????????
26
4functional_3/sequential/functional_1/flatten/Reshape?
5functional_3/sequential/functional_1/dropout/IdentityIdentity=functional_3/sequential/functional_1/flatten/Reshape:output:0*
T0*(
_output_shapes
:??????????
27
5functional_3/sequential/functional_1/dropout/Identity?
Afunctional_3/sequential/functional_1/layer1/MatMul/ReadVariableOpReadVariableOpJfunctional_3_sequential_functional_1_layer1_matmul_readvariableop_resource* 
_output_shapes
:
?
?*
dtype02C
Afunctional_3/sequential/functional_1/layer1/MatMul/ReadVariableOp?
2functional_3/sequential/functional_1/layer1/MatMulMatMul>functional_3/sequential/functional_1/dropout/Identity:output:0Ifunctional_3/sequential/functional_1/layer1/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????24
2functional_3/sequential/functional_1/layer1/MatMul?
Bfunctional_3/sequential/functional_1/layer1/BiasAdd/ReadVariableOpReadVariableOpKfunctional_3_sequential_functional_1_layer1_biasadd_readvariableop_resource*
_output_shapes	
:?*
dtype02D
Bfunctional_3/sequential/functional_1/layer1/BiasAdd/ReadVariableOp?
3functional_3/sequential/functional_1/layer1/BiasAddBiasAdd<functional_3/sequential/functional_1/layer1/MatMul:product:0Jfunctional_3/sequential/functional_1/layer1/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????25
3functional_3/sequential/functional_1/layer1/BiasAdd?
/functional_3/sequential/functional_1/layer1/EluElu<functional_3/sequential/functional_1/layer1/BiasAdd:output:0*
T0*(
_output_shapes
:??????????21
/functional_3/sequential/functional_1/layer1/Elu?
3functional_3/sequential/dense/MatMul/ReadVariableOpReadVariableOp<functional_3_sequential_dense_matmul_readvariableop_resource*
_output_shapes
:	?d*
dtype025
3functional_3/sequential/dense/MatMul/ReadVariableOp?
$functional_3/sequential/dense/MatMulMatMul=functional_3/sequential/functional_1/layer1/Elu:activations:0;functional_3/sequential/dense/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????d2&
$functional_3/sequential/dense/MatMul?
4functional_3/sequential/dense/BiasAdd/ReadVariableOpReadVariableOp=functional_3_sequential_dense_biasadd_readvariableop_resource*
_output_shapes
:d*
dtype026
4functional_3/sequential/dense/BiasAdd/ReadVariableOp?
%functional_3/sequential/dense/BiasAddBiasAdd.functional_3/sequential/dense/MatMul:product:0<functional_3/sequential/dense/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????d2'
%functional_3/sequential/dense/BiasAdd?
!functional_3/sequential/dense/EluElu.functional_3/sequential/dense/BiasAdd:output:0*
T0*'
_output_shapes
:?????????d2#
!functional_3/sequential/dense/Elu?
5functional_3/sequential/dense_1/MatMul/ReadVariableOpReadVariableOp>functional_3_sequential_dense_1_matmul_readvariableop_resource*
_output_shapes

:d
*
dtype027
5functional_3/sequential/dense_1/MatMul/ReadVariableOp?
&functional_3/sequential/dense_1/MatMulMatMul/functional_3/sequential/dense/Elu:activations:0=functional_3/sequential/dense_1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????
2(
&functional_3/sequential/dense_1/MatMul?
6functional_3/sequential/dense_1/BiasAdd/ReadVariableOpReadVariableOp?functional_3_sequential_dense_1_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype028
6functional_3/sequential/dense_1/BiasAdd/ReadVariableOp?
'functional_3/sequential/dense_1/BiasAddBiasAdd0functional_3/sequential/dense_1/MatMul:product:0>functional_3/sequential/dense_1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????
2)
'functional_3/sequential/dense_1/BiasAdd?
#functional_3/sequential/dense_1/EluElu0functional_3/sequential/dense_1/BiasAdd:output:0*
T0*'
_output_shapes
:?????????
2%
#functional_3/sequential/dense_1/Elu?
Cfunctional_3/sequential/functional_1/conv3d/Conv3D_1/ReadVariableOpReadVariableOpJfunctional_3_sequential_functional_1_conv3d_conv3d_readvariableop_resource*+
_output_shapes
:?*
dtype02E
Cfunctional_3/sequential/functional_1/conv3d/Conv3D_1/ReadVariableOp?
4functional_3/sequential/functional_1/conv3d/Conv3D_1Conv3Dinput_2Kfunctional_3/sequential/functional_1/conv3d/Conv3D_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
26
4functional_3/sequential/functional_1/conv3d/Conv3D_1?
Dfunctional_3/sequential/functional_1/conv3d/BiasAdd_1/ReadVariableOpReadVariableOpKfunctional_3_sequential_functional_1_conv3d_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02F
Dfunctional_3/sequential/functional_1/conv3d/BiasAdd_1/ReadVariableOp?
5functional_3/sequential/functional_1/conv3d/BiasAdd_1BiasAdd=functional_3/sequential/functional_1/conv3d/Conv3D_1:output:0Lfunctional_3/sequential/functional_1/conv3d/BiasAdd_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????27
5functional_3/sequential/functional_1/conv3d/BiasAdd_1?
Sfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/ReadVariableOpReadVariableOpZfunctional_3_sequential_functional_1_batch_normalization_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02U
Sfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/ReadVariableOp?
Jfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2L
Jfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/add/y?
Hfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/addAddV2[functional_3/sequential/functional_1/batch_normalization/batchnorm_1/ReadVariableOp:value:0Sfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/add/y:output:0*
T0*
_output_shapes
:2J
Hfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/add?
Jfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/RsqrtRsqrtLfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/add:z:0*
T0*
_output_shapes
:2L
Jfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/Rsqrt?
Wfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/mul/ReadVariableOpReadVariableOp^functional_3_sequential_functional_1_batch_normalization_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02Y
Wfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/mul/ReadVariableOp?
Hfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/mulMulNfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/Rsqrt:y:0_functional_3/sequential/functional_1/batch_normalization/batchnorm_1/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2J
Hfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/mul?
Jfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/mul_1Mul>functional_3/sequential/functional_1/conv3d/BiasAdd_1:output:0Lfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/mul:z:0*
T0*3
_output_shapes!
:?????????2L
Jfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/mul_1?
Ufunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/ReadVariableOp_1ReadVariableOp\functional_3_sequential_functional_1_batch_normalization_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02W
Ufunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/ReadVariableOp_1?
Jfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/mul_2Mul]functional_3/sequential/functional_1/batch_normalization/batchnorm_1/ReadVariableOp_1:value:0Lfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/mul:z:0*
T0*
_output_shapes
:2L
Jfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/mul_2?
Ufunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/ReadVariableOp_2ReadVariableOp\functional_3_sequential_functional_1_batch_normalization_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02W
Ufunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/ReadVariableOp_2?
Hfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/subSub]functional_3/sequential/functional_1/batch_normalization/batchnorm_1/ReadVariableOp_2:value:0Nfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/mul_2:z:0*
T0*
_output_shapes
:2J
Hfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/sub?
Jfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/add_1AddV2Nfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/mul_1:z:0Lfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/sub:z:0*
T0*3
_output_shapes!
:?????????2L
Jfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/add_1?
Efunctional_3/sequential/functional_1/conv3d_1/Conv3D_1/ReadVariableOpReadVariableOpLfunctional_3_sequential_functional_1_conv3d_1_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02G
Efunctional_3/sequential/functional_1/conv3d_1/Conv3D_1/ReadVariableOp?
6functional_3/sequential/functional_1/conv3d_1/Conv3D_1Conv3DNfunctional_3/sequential/functional_1/batch_normalization/batchnorm_1/add_1:z:0Mfunctional_3/sequential/functional_1/conv3d_1/Conv3D_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
28
6functional_3/sequential/functional_1/conv3d_1/Conv3D_1?
Ffunctional_3/sequential/functional_1/conv3d_1/BiasAdd_1/ReadVariableOpReadVariableOpMfunctional_3_sequential_functional_1_conv3d_1_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02H
Ffunctional_3/sequential/functional_1/conv3d_1/BiasAdd_1/ReadVariableOp?
7functional_3/sequential/functional_1/conv3d_1/BiasAdd_1BiasAdd?functional_3/sequential/functional_1/conv3d_1/Conv3D_1:output:0Nfunctional_3/sequential/functional_1/conv3d_1/BiasAdd_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????29
7functional_3/sequential/functional_1/conv3d_1/BiasAdd_1?
3functional_3/sequential/functional_1/conv3d_1/Elu_1Elu@functional_3/sequential/functional_1/conv3d_1/BiasAdd_1:output:0*
T0*3
_output_shapes!
:?????????25
3functional_3/sequential/functional_1/conv3d_1/Elu_1?
Ufunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/ReadVariableOpReadVariableOp\functional_3_sequential_functional_1_batch_normalization_1_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02W
Ufunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/ReadVariableOp?
Lfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2N
Lfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/add/y?
Jfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/addAddV2]functional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/ReadVariableOp:value:0Ufunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/add/y:output:0*
T0*
_output_shapes
:2L
Jfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/add?
Lfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/RsqrtRsqrtNfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/add:z:0*
T0*
_output_shapes
:2N
Lfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/Rsqrt?
Yfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/mul/ReadVariableOpReadVariableOp`functional_3_sequential_functional_1_batch_normalization_1_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02[
Yfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/mul/ReadVariableOp?
Jfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/mulMulPfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/Rsqrt:y:0afunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2L
Jfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/mul?
Lfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/mul_1MulAfunctional_3/sequential/functional_1/conv3d_1/Elu_1:activations:0Nfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/mul:z:0*
T0*3
_output_shapes!
:?????????2N
Lfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/mul_1?
Wfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/ReadVariableOp_1ReadVariableOp^functional_3_sequential_functional_1_batch_normalization_1_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02Y
Wfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/ReadVariableOp_1?
Lfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/mul_2Mul_functional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/ReadVariableOp_1:value:0Nfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/mul:z:0*
T0*
_output_shapes
:2N
Lfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/mul_2?
Wfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/ReadVariableOp_2ReadVariableOp^functional_3_sequential_functional_1_batch_normalization_1_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02Y
Wfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/ReadVariableOp_2?
Jfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/subSub_functional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/ReadVariableOp_2:value:0Pfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/mul_2:z:0*
T0*
_output_shapes
:2L
Jfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/sub?
Lfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/add_1AddV2Pfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/mul_1:z:0Nfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/sub:z:0*
T0*3
_output_shapes!
:?????????2N
Lfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/add_1?
Efunctional_3/sequential/functional_1/conv3d_2/Conv3D_1/ReadVariableOpReadVariableOpLfunctional_3_sequential_functional_1_conv3d_2_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02G
Efunctional_3/sequential/functional_1/conv3d_2/Conv3D_1/ReadVariableOp?
6functional_3/sequential/functional_1/conv3d_2/Conv3D_1Conv3DPfunctional_3/sequential/functional_1/batch_normalization_1/batchnorm_1/add_1:z:0Mfunctional_3/sequential/functional_1/conv3d_2/Conv3D_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
28
6functional_3/sequential/functional_1/conv3d_2/Conv3D_1?
Ffunctional_3/sequential/functional_1/conv3d_2/BiasAdd_1/ReadVariableOpReadVariableOpMfunctional_3_sequential_functional_1_conv3d_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02H
Ffunctional_3/sequential/functional_1/conv3d_2/BiasAdd_1/ReadVariableOp?
7functional_3/sequential/functional_1/conv3d_2/BiasAdd_1BiasAdd?functional_3/sequential/functional_1/conv3d_2/Conv3D_1:output:0Nfunctional_3/sequential/functional_1/conv3d_2/BiasAdd_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????29
7functional_3/sequential/functional_1/conv3d_2/BiasAdd_1?
3functional_3/sequential/functional_1/conv3d_2/Elu_1Elu@functional_3/sequential/functional_1/conv3d_2/BiasAdd_1:output:0*
T0*3
_output_shapes!
:?????????25
3functional_3/sequential/functional_1/conv3d_2/Elu_1?
Ufunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/ReadVariableOpReadVariableOp\functional_3_sequential_functional_1_batch_normalization_2_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02W
Ufunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/ReadVariableOp?
Lfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2N
Lfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/add/y?
Jfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/addAddV2]functional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/ReadVariableOp:value:0Ufunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/add/y:output:0*
T0*
_output_shapes
:2L
Jfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/add?
Lfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/RsqrtRsqrtNfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/add:z:0*
T0*
_output_shapes
:2N
Lfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/Rsqrt?
Yfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/mul/ReadVariableOpReadVariableOp`functional_3_sequential_functional_1_batch_normalization_2_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02[
Yfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/mul/ReadVariableOp?
Jfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/mulMulPfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/Rsqrt:y:0afunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2L
Jfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/mul?
Lfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/mul_1MulAfunctional_3/sequential/functional_1/conv3d_2/Elu_1:activations:0Nfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/mul:z:0*
T0*3
_output_shapes!
:?????????2N
Lfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/mul_1?
Wfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/ReadVariableOp_1ReadVariableOp^functional_3_sequential_functional_1_batch_normalization_2_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02Y
Wfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/ReadVariableOp_1?
Lfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/mul_2Mul_functional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/ReadVariableOp_1:value:0Nfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/mul:z:0*
T0*
_output_shapes
:2N
Lfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/mul_2?
Wfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/ReadVariableOp_2ReadVariableOp^functional_3_sequential_functional_1_batch_normalization_2_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02Y
Wfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/ReadVariableOp_2?
Jfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/subSub_functional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/ReadVariableOp_2:value:0Pfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/mul_2:z:0*
T0*
_output_shapes
:2L
Jfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/sub?
Lfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/add_1AddV2Pfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/mul_1:z:0Nfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/sub:z:0*
T0*3
_output_shapes!
:?????????2N
Lfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/add_1?
Efunctional_3/sequential/functional_1/conv3d_3/Conv3D_1/ReadVariableOpReadVariableOpLfunctional_3_sequential_functional_1_conv3d_3_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02G
Efunctional_3/sequential/functional_1/conv3d_3/Conv3D_1/ReadVariableOp?
6functional_3/sequential/functional_1/conv3d_3/Conv3D_1Conv3DPfunctional_3/sequential/functional_1/batch_normalization_2/batchnorm_1/add_1:z:0Mfunctional_3/sequential/functional_1/conv3d_3/Conv3D_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
28
6functional_3/sequential/functional_1/conv3d_3/Conv3D_1?
Ffunctional_3/sequential/functional_1/conv3d_3/BiasAdd_1/ReadVariableOpReadVariableOpMfunctional_3_sequential_functional_1_conv3d_3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02H
Ffunctional_3/sequential/functional_1/conv3d_3/BiasAdd_1/ReadVariableOp?
7functional_3/sequential/functional_1/conv3d_3/BiasAdd_1BiasAdd?functional_3/sequential/functional_1/conv3d_3/Conv3D_1:output:0Nfunctional_3/sequential/functional_1/conv3d_3/BiasAdd_1/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????29
7functional_3/sequential/functional_1/conv3d_3/BiasAdd_1?
3functional_3/sequential/functional_1/conv3d_3/Elu_1Elu@functional_3/sequential/functional_1/conv3d_3/BiasAdd_1:output:0*
T0*3
_output_shapes!
:?????????25
3functional_3/sequential/functional_1/conv3d_3/Elu_1?
Ufunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/ReadVariableOpReadVariableOp\functional_3_sequential_functional_1_batch_normalization_3_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02W
Ufunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/ReadVariableOp?
Lfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2N
Lfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/add/y?
Jfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/addAddV2]functional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/ReadVariableOp:value:0Ufunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/add/y:output:0*
T0*
_output_shapes
:2L
Jfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/add?
Lfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/RsqrtRsqrtNfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/add:z:0*
T0*
_output_shapes
:2N
Lfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/Rsqrt?
Yfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/mul/ReadVariableOpReadVariableOp`functional_3_sequential_functional_1_batch_normalization_3_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02[
Yfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/mul/ReadVariableOp?
Jfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/mulMulPfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/Rsqrt:y:0afunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2L
Jfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/mul?
Lfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/mul_1MulAfunctional_3/sequential/functional_1/conv3d_3/Elu_1:activations:0Nfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/mul:z:0*
T0*3
_output_shapes!
:?????????2N
Lfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/mul_1?
Wfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/ReadVariableOp_1ReadVariableOp^functional_3_sequential_functional_1_batch_normalization_3_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02Y
Wfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/ReadVariableOp_1?
Lfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/mul_2Mul_functional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/ReadVariableOp_1:value:0Nfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/mul:z:0*
T0*
_output_shapes
:2N
Lfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/mul_2?
Wfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/ReadVariableOp_2ReadVariableOp^functional_3_sequential_functional_1_batch_normalization_3_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02Y
Wfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/ReadVariableOp_2?
Jfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/subSub_functional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/ReadVariableOp_2:value:0Pfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/mul_2:z:0*
T0*
_output_shapes
:2L
Jfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/sub?
Lfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/add_1AddV2Pfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/mul_1:z:0Nfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/sub:z:0*
T0*3
_output_shapes!
:?????????2N
Lfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/add_1?
Bfunctional_3/sequential/functional_1/average_pooling3d/AvgPool3D_1	AvgPool3DPfunctional_3/sequential/functional_1/batch_normalization_3/batchnorm_1/add_1:z:0*
T0*3
_output_shapes!
:?????????*
ksize	
*
paddingVALID*
strides	
2D
Bfunctional_3/sequential/functional_1/average_pooling3d/AvgPool3D_1?
4functional_3/sequential/functional_1/flatten/Const_1Const*
_output_shapes
:*
dtype0*
valueB"????   26
4functional_3/sequential/functional_1/flatten/Const_1?
6functional_3/sequential/functional_1/flatten/Reshape_1ReshapeKfunctional_3/sequential/functional_1/average_pooling3d/AvgPool3D_1:output:0=functional_3/sequential/functional_1/flatten/Const_1:output:0*
T0*(
_output_shapes
:??????????
28
6functional_3/sequential/functional_1/flatten/Reshape_1?
7functional_3/sequential/functional_1/dropout/Identity_1Identity?functional_3/sequential/functional_1/flatten/Reshape_1:output:0*
T0*(
_output_shapes
:??????????
29
7functional_3/sequential/functional_1/dropout/Identity_1?
Cfunctional_3/sequential/functional_1/layer1/MatMul_1/ReadVariableOpReadVariableOpJfunctional_3_sequential_functional_1_layer1_matmul_readvariableop_resource* 
_output_shapes
:
?
?*
dtype02E
Cfunctional_3/sequential/functional_1/layer1/MatMul_1/ReadVariableOp?
4functional_3/sequential/functional_1/layer1/MatMul_1MatMul@functional_3/sequential/functional_1/dropout/Identity_1:output:0Kfunctional_3/sequential/functional_1/layer1/MatMul_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????26
4functional_3/sequential/functional_1/layer1/MatMul_1?
Dfunctional_3/sequential/functional_1/layer1/BiasAdd_1/ReadVariableOpReadVariableOpKfunctional_3_sequential_functional_1_layer1_biasadd_readvariableop_resource*
_output_shapes	
:?*
dtype02F
Dfunctional_3/sequential/functional_1/layer1/BiasAdd_1/ReadVariableOp?
5functional_3/sequential/functional_1/layer1/BiasAdd_1BiasAdd>functional_3/sequential/functional_1/layer1/MatMul_1:product:0Lfunctional_3/sequential/functional_1/layer1/BiasAdd_1/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????27
5functional_3/sequential/functional_1/layer1/BiasAdd_1?
1functional_3/sequential/functional_1/layer1/Elu_1Elu>functional_3/sequential/functional_1/layer1/BiasAdd_1:output:0*
T0*(
_output_shapes
:??????????23
1functional_3/sequential/functional_1/layer1/Elu_1?
5functional_3/sequential/dense/MatMul_1/ReadVariableOpReadVariableOp<functional_3_sequential_dense_matmul_readvariableop_resource*
_output_shapes
:	?d*
dtype027
5functional_3/sequential/dense/MatMul_1/ReadVariableOp?
&functional_3/sequential/dense/MatMul_1MatMul?functional_3/sequential/functional_1/layer1/Elu_1:activations:0=functional_3/sequential/dense/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????d2(
&functional_3/sequential/dense/MatMul_1?
6functional_3/sequential/dense/BiasAdd_1/ReadVariableOpReadVariableOp=functional_3_sequential_dense_biasadd_readvariableop_resource*
_output_shapes
:d*
dtype028
6functional_3/sequential/dense/BiasAdd_1/ReadVariableOp?
'functional_3/sequential/dense/BiasAdd_1BiasAdd0functional_3/sequential/dense/MatMul_1:product:0>functional_3/sequential/dense/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????d2)
'functional_3/sequential/dense/BiasAdd_1?
#functional_3/sequential/dense/Elu_1Elu0functional_3/sequential/dense/BiasAdd_1:output:0*
T0*'
_output_shapes
:?????????d2%
#functional_3/sequential/dense/Elu_1?
7functional_3/sequential/dense_1/MatMul_1/ReadVariableOpReadVariableOp>functional_3_sequential_dense_1_matmul_readvariableop_resource*
_output_shapes

:d
*
dtype029
7functional_3/sequential/dense_1/MatMul_1/ReadVariableOp?
(functional_3/sequential/dense_1/MatMul_1MatMul1functional_3/sequential/dense/Elu_1:activations:0?functional_3/sequential/dense_1/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????
2*
(functional_3/sequential/dense_1/MatMul_1?
8functional_3/sequential/dense_1/BiasAdd_1/ReadVariableOpReadVariableOp?functional_3_sequential_dense_1_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02:
8functional_3/sequential/dense_1/BiasAdd_1/ReadVariableOp?
)functional_3/sequential/dense_1/BiasAdd_1BiasAdd2functional_3/sequential/dense_1/MatMul_1:product:0@functional_3/sequential/dense_1/BiasAdd_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????
2+
)functional_3/sequential/dense_1/BiasAdd_1?
%functional_3/sequential/dense_1/Elu_1Elu2functional_3/sequential/dense_1/BiasAdd_1:output:0*
T0*'
_output_shapes
:?????????
2'
%functional_3/sequential/dense_1/Elu_1?
functional_3/lambda/subSub1functional_3/sequential/dense_1/Elu:activations:03functional_3/sequential/dense_1/Elu_1:activations:0*
T0*'
_output_shapes
:?????????
2
functional_3/lambda/sub?
functional_3/lambda/AbsAbsfunctional_3/lambda/sub:z:0*
T0*'
_output_shapes
:?????????
2
functional_3/lambda/Abs?
$functional_3/concatenate/concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2&
$functional_3/concatenate/concat/axis?
functional_3/concatenate/concatConcatV2functional_3/lambda/Abs:y:0input_3-functional_3/concatenate/concat/axis:output:0*
N*
T0*'
_output_shapes
:?????????2!
functional_3/concatenate/concat?
*functional_3/dense_2/MatMul/ReadVariableOpReadVariableOp3functional_3_dense_2_matmul_readvariableop_resource*
_output_shapes

:*
dtype02,
*functional_3/dense_2/MatMul/ReadVariableOp?
functional_3/dense_2/MatMulMatMul(functional_3/concatenate/concat:output:02functional_3/dense_2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
functional_3/dense_2/MatMul?
+functional_3/dense_2/BiasAdd/ReadVariableOpReadVariableOp4functional_3_dense_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02-
+functional_3/dense_2/BiasAdd/ReadVariableOp?
functional_3/dense_2/BiasAddBiasAdd%functional_3/dense_2/MatMul:product:03functional_3/dense_2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
functional_3/dense_2/BiasAddy
IdentityIdentity%functional_3/dense_2/BiasAdd:output:0*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????:??????????:?????????
:::::::::::::::::::::::::::::::::] Y
4
_output_shapes"
 :??????????
!
_user_specified_name	input_1:]Y
4
_output_shapes"
 :??????????
!
_user_specified_name	input_2:PL
'
_output_shapes
:?????????

!
_user_specified_name	input_3
??
?
G__inference_sequential_layer_call_and_return_conditional_losses_1226239

inputs6
2functional_1_conv3d_conv3d_readvariableop_resource7
3functional_1_conv3d_biasadd_readvariableop_resourceF
Bfunctional_1_batch_normalization_batchnorm_readvariableop_resourceJ
Ffunctional_1_batch_normalization_batchnorm_mul_readvariableop_resourceH
Dfunctional_1_batch_normalization_batchnorm_readvariableop_1_resourceH
Dfunctional_1_batch_normalization_batchnorm_readvariableop_2_resource8
4functional_1_conv3d_1_conv3d_readvariableop_resource9
5functional_1_conv3d_1_biasadd_readvariableop_resourceH
Dfunctional_1_batch_normalization_1_batchnorm_readvariableop_resourceL
Hfunctional_1_batch_normalization_1_batchnorm_mul_readvariableop_resourceJ
Ffunctional_1_batch_normalization_1_batchnorm_readvariableop_1_resourceJ
Ffunctional_1_batch_normalization_1_batchnorm_readvariableop_2_resource8
4functional_1_conv3d_2_conv3d_readvariableop_resource9
5functional_1_conv3d_2_biasadd_readvariableop_resourceH
Dfunctional_1_batch_normalization_2_batchnorm_readvariableop_resourceL
Hfunctional_1_batch_normalization_2_batchnorm_mul_readvariableop_resourceJ
Ffunctional_1_batch_normalization_2_batchnorm_readvariableop_1_resourceJ
Ffunctional_1_batch_normalization_2_batchnorm_readvariableop_2_resource8
4functional_1_conv3d_3_conv3d_readvariableop_resource9
5functional_1_conv3d_3_biasadd_readvariableop_resourceH
Dfunctional_1_batch_normalization_3_batchnorm_readvariableop_resourceL
Hfunctional_1_batch_normalization_3_batchnorm_mul_readvariableop_resourceJ
Ffunctional_1_batch_normalization_3_batchnorm_readvariableop_1_resourceJ
Ffunctional_1_batch_normalization_3_batchnorm_readvariableop_2_resource6
2functional_1_layer1_matmul_readvariableop_resource7
3functional_1_layer1_biasadd_readvariableop_resource(
$dense_matmul_readvariableop_resource)
%dense_biasadd_readvariableop_resource*
&dense_1_matmul_readvariableop_resource+
'dense_1_biasadd_readvariableop_resource
identity??
)functional_1/conv3d/Conv3D/ReadVariableOpReadVariableOp2functional_1_conv3d_conv3d_readvariableop_resource*+
_output_shapes
:?*
dtype02+
)functional_1/conv3d/Conv3D/ReadVariableOp?
functional_1/conv3d/Conv3DConv3Dinputs1functional_1/conv3d/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
functional_1/conv3d/Conv3D?
*functional_1/conv3d/BiasAdd/ReadVariableOpReadVariableOp3functional_1_conv3d_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02,
*functional_1/conv3d/BiasAdd/ReadVariableOp?
functional_1/conv3d/BiasAddBiasAdd#functional_1/conv3d/Conv3D:output:02functional_1/conv3d/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2
functional_1/conv3d/BiasAdd?
9functional_1/batch_normalization/batchnorm/ReadVariableOpReadVariableOpBfunctional_1_batch_normalization_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02;
9functional_1/batch_normalization/batchnorm/ReadVariableOp?
0functional_1/batch_normalization/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:22
0functional_1/batch_normalization/batchnorm/add/y?
.functional_1/batch_normalization/batchnorm/addAddV2Afunctional_1/batch_normalization/batchnorm/ReadVariableOp:value:09functional_1/batch_normalization/batchnorm/add/y:output:0*
T0*
_output_shapes
:20
.functional_1/batch_normalization/batchnorm/add?
0functional_1/batch_normalization/batchnorm/RsqrtRsqrt2functional_1/batch_normalization/batchnorm/add:z:0*
T0*
_output_shapes
:22
0functional_1/batch_normalization/batchnorm/Rsqrt?
=functional_1/batch_normalization/batchnorm/mul/ReadVariableOpReadVariableOpFfunctional_1_batch_normalization_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02?
=functional_1/batch_normalization/batchnorm/mul/ReadVariableOp?
.functional_1/batch_normalization/batchnorm/mulMul4functional_1/batch_normalization/batchnorm/Rsqrt:y:0Efunctional_1/batch_normalization/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:20
.functional_1/batch_normalization/batchnorm/mul?
0functional_1/batch_normalization/batchnorm/mul_1Mul$functional_1/conv3d/BiasAdd:output:02functional_1/batch_normalization/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????22
0functional_1/batch_normalization/batchnorm/mul_1?
;functional_1/batch_normalization/batchnorm/ReadVariableOp_1ReadVariableOpDfunctional_1_batch_normalization_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02=
;functional_1/batch_normalization/batchnorm/ReadVariableOp_1?
0functional_1/batch_normalization/batchnorm/mul_2MulCfunctional_1/batch_normalization/batchnorm/ReadVariableOp_1:value:02functional_1/batch_normalization/batchnorm/mul:z:0*
T0*
_output_shapes
:22
0functional_1/batch_normalization/batchnorm/mul_2?
;functional_1/batch_normalization/batchnorm/ReadVariableOp_2ReadVariableOpDfunctional_1_batch_normalization_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02=
;functional_1/batch_normalization/batchnorm/ReadVariableOp_2?
.functional_1/batch_normalization/batchnorm/subSubCfunctional_1/batch_normalization/batchnorm/ReadVariableOp_2:value:04functional_1/batch_normalization/batchnorm/mul_2:z:0*
T0*
_output_shapes
:20
.functional_1/batch_normalization/batchnorm/sub?
0functional_1/batch_normalization/batchnorm/add_1AddV24functional_1/batch_normalization/batchnorm/mul_1:z:02functional_1/batch_normalization/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????22
0functional_1/batch_normalization/batchnorm/add_1?
+functional_1/conv3d_1/Conv3D/ReadVariableOpReadVariableOp4functional_1_conv3d_1_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02-
+functional_1/conv3d_1/Conv3D/ReadVariableOp?
functional_1/conv3d_1/Conv3DConv3D4functional_1/batch_normalization/batchnorm/add_1:z:03functional_1/conv3d_1/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
functional_1/conv3d_1/Conv3D?
,functional_1/conv3d_1/BiasAdd/ReadVariableOpReadVariableOp5functional_1_conv3d_1_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02.
,functional_1/conv3d_1/BiasAdd/ReadVariableOp?
functional_1/conv3d_1/BiasAddBiasAdd%functional_1/conv3d_1/Conv3D:output:04functional_1/conv3d_1/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2
functional_1/conv3d_1/BiasAdd?
functional_1/conv3d_1/EluElu&functional_1/conv3d_1/BiasAdd:output:0*
T0*3
_output_shapes!
:?????????2
functional_1/conv3d_1/Elu?
;functional_1/batch_normalization_1/batchnorm/ReadVariableOpReadVariableOpDfunctional_1_batch_normalization_1_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02=
;functional_1/batch_normalization_1/batchnorm/ReadVariableOp?
2functional_1/batch_normalization_1/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:24
2functional_1/batch_normalization_1/batchnorm/add/y?
0functional_1/batch_normalization_1/batchnorm/addAddV2Cfunctional_1/batch_normalization_1/batchnorm/ReadVariableOp:value:0;functional_1/batch_normalization_1/batchnorm/add/y:output:0*
T0*
_output_shapes
:22
0functional_1/batch_normalization_1/batchnorm/add?
2functional_1/batch_normalization_1/batchnorm/RsqrtRsqrt4functional_1/batch_normalization_1/batchnorm/add:z:0*
T0*
_output_shapes
:24
2functional_1/batch_normalization_1/batchnorm/Rsqrt?
?functional_1/batch_normalization_1/batchnorm/mul/ReadVariableOpReadVariableOpHfunctional_1_batch_normalization_1_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02A
?functional_1/batch_normalization_1/batchnorm/mul/ReadVariableOp?
0functional_1/batch_normalization_1/batchnorm/mulMul6functional_1/batch_normalization_1/batchnorm/Rsqrt:y:0Gfunctional_1/batch_normalization_1/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:22
0functional_1/batch_normalization_1/batchnorm/mul?
2functional_1/batch_normalization_1/batchnorm/mul_1Mul'functional_1/conv3d_1/Elu:activations:04functional_1/batch_normalization_1/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????24
2functional_1/batch_normalization_1/batchnorm/mul_1?
=functional_1/batch_normalization_1/batchnorm/ReadVariableOp_1ReadVariableOpFfunctional_1_batch_normalization_1_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02?
=functional_1/batch_normalization_1/batchnorm/ReadVariableOp_1?
2functional_1/batch_normalization_1/batchnorm/mul_2MulEfunctional_1/batch_normalization_1/batchnorm/ReadVariableOp_1:value:04functional_1/batch_normalization_1/batchnorm/mul:z:0*
T0*
_output_shapes
:24
2functional_1/batch_normalization_1/batchnorm/mul_2?
=functional_1/batch_normalization_1/batchnorm/ReadVariableOp_2ReadVariableOpFfunctional_1_batch_normalization_1_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02?
=functional_1/batch_normalization_1/batchnorm/ReadVariableOp_2?
0functional_1/batch_normalization_1/batchnorm/subSubEfunctional_1/batch_normalization_1/batchnorm/ReadVariableOp_2:value:06functional_1/batch_normalization_1/batchnorm/mul_2:z:0*
T0*
_output_shapes
:22
0functional_1/batch_normalization_1/batchnorm/sub?
2functional_1/batch_normalization_1/batchnorm/add_1AddV26functional_1/batch_normalization_1/batchnorm/mul_1:z:04functional_1/batch_normalization_1/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????24
2functional_1/batch_normalization_1/batchnorm/add_1?
+functional_1/conv3d_2/Conv3D/ReadVariableOpReadVariableOp4functional_1_conv3d_2_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02-
+functional_1/conv3d_2/Conv3D/ReadVariableOp?
functional_1/conv3d_2/Conv3DConv3D6functional_1/batch_normalization_1/batchnorm/add_1:z:03functional_1/conv3d_2/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
functional_1/conv3d_2/Conv3D?
,functional_1/conv3d_2/BiasAdd/ReadVariableOpReadVariableOp5functional_1_conv3d_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02.
,functional_1/conv3d_2/BiasAdd/ReadVariableOp?
functional_1/conv3d_2/BiasAddBiasAdd%functional_1/conv3d_2/Conv3D:output:04functional_1/conv3d_2/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2
functional_1/conv3d_2/BiasAdd?
functional_1/conv3d_2/EluElu&functional_1/conv3d_2/BiasAdd:output:0*
T0*3
_output_shapes!
:?????????2
functional_1/conv3d_2/Elu?
;functional_1/batch_normalization_2/batchnorm/ReadVariableOpReadVariableOpDfunctional_1_batch_normalization_2_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02=
;functional_1/batch_normalization_2/batchnorm/ReadVariableOp?
2functional_1/batch_normalization_2/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:24
2functional_1/batch_normalization_2/batchnorm/add/y?
0functional_1/batch_normalization_2/batchnorm/addAddV2Cfunctional_1/batch_normalization_2/batchnorm/ReadVariableOp:value:0;functional_1/batch_normalization_2/batchnorm/add/y:output:0*
T0*
_output_shapes
:22
0functional_1/batch_normalization_2/batchnorm/add?
2functional_1/batch_normalization_2/batchnorm/RsqrtRsqrt4functional_1/batch_normalization_2/batchnorm/add:z:0*
T0*
_output_shapes
:24
2functional_1/batch_normalization_2/batchnorm/Rsqrt?
?functional_1/batch_normalization_2/batchnorm/mul/ReadVariableOpReadVariableOpHfunctional_1_batch_normalization_2_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02A
?functional_1/batch_normalization_2/batchnorm/mul/ReadVariableOp?
0functional_1/batch_normalization_2/batchnorm/mulMul6functional_1/batch_normalization_2/batchnorm/Rsqrt:y:0Gfunctional_1/batch_normalization_2/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:22
0functional_1/batch_normalization_2/batchnorm/mul?
2functional_1/batch_normalization_2/batchnorm/mul_1Mul'functional_1/conv3d_2/Elu:activations:04functional_1/batch_normalization_2/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????24
2functional_1/batch_normalization_2/batchnorm/mul_1?
=functional_1/batch_normalization_2/batchnorm/ReadVariableOp_1ReadVariableOpFfunctional_1_batch_normalization_2_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02?
=functional_1/batch_normalization_2/batchnorm/ReadVariableOp_1?
2functional_1/batch_normalization_2/batchnorm/mul_2MulEfunctional_1/batch_normalization_2/batchnorm/ReadVariableOp_1:value:04functional_1/batch_normalization_2/batchnorm/mul:z:0*
T0*
_output_shapes
:24
2functional_1/batch_normalization_2/batchnorm/mul_2?
=functional_1/batch_normalization_2/batchnorm/ReadVariableOp_2ReadVariableOpFfunctional_1_batch_normalization_2_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02?
=functional_1/batch_normalization_2/batchnorm/ReadVariableOp_2?
0functional_1/batch_normalization_2/batchnorm/subSubEfunctional_1/batch_normalization_2/batchnorm/ReadVariableOp_2:value:06functional_1/batch_normalization_2/batchnorm/mul_2:z:0*
T0*
_output_shapes
:22
0functional_1/batch_normalization_2/batchnorm/sub?
2functional_1/batch_normalization_2/batchnorm/add_1AddV26functional_1/batch_normalization_2/batchnorm/mul_1:z:04functional_1/batch_normalization_2/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????24
2functional_1/batch_normalization_2/batchnorm/add_1?
+functional_1/conv3d_3/Conv3D/ReadVariableOpReadVariableOp4functional_1_conv3d_3_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02-
+functional_1/conv3d_3/Conv3D/ReadVariableOp?
functional_1/conv3d_3/Conv3DConv3D6functional_1/batch_normalization_2/batchnorm/add_1:z:03functional_1/conv3d_3/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
functional_1/conv3d_3/Conv3D?
,functional_1/conv3d_3/BiasAdd/ReadVariableOpReadVariableOp5functional_1_conv3d_3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02.
,functional_1/conv3d_3/BiasAdd/ReadVariableOp?
functional_1/conv3d_3/BiasAddBiasAdd%functional_1/conv3d_3/Conv3D:output:04functional_1/conv3d_3/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2
functional_1/conv3d_3/BiasAdd?
functional_1/conv3d_3/EluElu&functional_1/conv3d_3/BiasAdd:output:0*
T0*3
_output_shapes!
:?????????2
functional_1/conv3d_3/Elu?
;functional_1/batch_normalization_3/batchnorm/ReadVariableOpReadVariableOpDfunctional_1_batch_normalization_3_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02=
;functional_1/batch_normalization_3/batchnorm/ReadVariableOp?
2functional_1/batch_normalization_3/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:24
2functional_1/batch_normalization_3/batchnorm/add/y?
0functional_1/batch_normalization_3/batchnorm/addAddV2Cfunctional_1/batch_normalization_3/batchnorm/ReadVariableOp:value:0;functional_1/batch_normalization_3/batchnorm/add/y:output:0*
T0*
_output_shapes
:22
0functional_1/batch_normalization_3/batchnorm/add?
2functional_1/batch_normalization_3/batchnorm/RsqrtRsqrt4functional_1/batch_normalization_3/batchnorm/add:z:0*
T0*
_output_shapes
:24
2functional_1/batch_normalization_3/batchnorm/Rsqrt?
?functional_1/batch_normalization_3/batchnorm/mul/ReadVariableOpReadVariableOpHfunctional_1_batch_normalization_3_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02A
?functional_1/batch_normalization_3/batchnorm/mul/ReadVariableOp?
0functional_1/batch_normalization_3/batchnorm/mulMul6functional_1/batch_normalization_3/batchnorm/Rsqrt:y:0Gfunctional_1/batch_normalization_3/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:22
0functional_1/batch_normalization_3/batchnorm/mul?
2functional_1/batch_normalization_3/batchnorm/mul_1Mul'functional_1/conv3d_3/Elu:activations:04functional_1/batch_normalization_3/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????24
2functional_1/batch_normalization_3/batchnorm/mul_1?
=functional_1/batch_normalization_3/batchnorm/ReadVariableOp_1ReadVariableOpFfunctional_1_batch_normalization_3_batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02?
=functional_1/batch_normalization_3/batchnorm/ReadVariableOp_1?
2functional_1/batch_normalization_3/batchnorm/mul_2MulEfunctional_1/batch_normalization_3/batchnorm/ReadVariableOp_1:value:04functional_1/batch_normalization_3/batchnorm/mul:z:0*
T0*
_output_shapes
:24
2functional_1/batch_normalization_3/batchnorm/mul_2?
=functional_1/batch_normalization_3/batchnorm/ReadVariableOp_2ReadVariableOpFfunctional_1_batch_normalization_3_batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02?
=functional_1/batch_normalization_3/batchnorm/ReadVariableOp_2?
0functional_1/batch_normalization_3/batchnorm/subSubEfunctional_1/batch_normalization_3/batchnorm/ReadVariableOp_2:value:06functional_1/batch_normalization_3/batchnorm/mul_2:z:0*
T0*
_output_shapes
:22
0functional_1/batch_normalization_3/batchnorm/sub?
2functional_1/batch_normalization_3/batchnorm/add_1AddV26functional_1/batch_normalization_3/batchnorm/mul_1:z:04functional_1/batch_normalization_3/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????24
2functional_1/batch_normalization_3/batchnorm/add_1?
(functional_1/average_pooling3d/AvgPool3D	AvgPool3D6functional_1/batch_normalization_3/batchnorm/add_1:z:0*
T0*3
_output_shapes!
:?????????*
ksize	
*
paddingVALID*
strides	
2*
(functional_1/average_pooling3d/AvgPool3D?
functional_1/flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"????   2
functional_1/flatten/Const?
functional_1/flatten/ReshapeReshape1functional_1/average_pooling3d/AvgPool3D:output:0#functional_1/flatten/Const:output:0*
T0*(
_output_shapes
:??????????
2
functional_1/flatten/Reshape?
functional_1/dropout/IdentityIdentity%functional_1/flatten/Reshape:output:0*
T0*(
_output_shapes
:??????????
2
functional_1/dropout/Identity?
)functional_1/layer1/MatMul/ReadVariableOpReadVariableOp2functional_1_layer1_matmul_readvariableop_resource* 
_output_shapes
:
?
?*
dtype02+
)functional_1/layer1/MatMul/ReadVariableOp?
functional_1/layer1/MatMulMatMul&functional_1/dropout/Identity:output:01functional_1/layer1/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2
functional_1/layer1/MatMul?
*functional_1/layer1/BiasAdd/ReadVariableOpReadVariableOp3functional_1_layer1_biasadd_readvariableop_resource*
_output_shapes	
:?*
dtype02,
*functional_1/layer1/BiasAdd/ReadVariableOp?
functional_1/layer1/BiasAddBiasAdd$functional_1/layer1/MatMul:product:02functional_1/layer1/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2
functional_1/layer1/BiasAdd?
functional_1/layer1/EluElu$functional_1/layer1/BiasAdd:output:0*
T0*(
_output_shapes
:??????????2
functional_1/layer1/Elu?
dense/MatMul/ReadVariableOpReadVariableOp$dense_matmul_readvariableop_resource*
_output_shapes
:	?d*
dtype02
dense/MatMul/ReadVariableOp?
dense/MatMulMatMul%functional_1/layer1/Elu:activations:0#dense/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????d2
dense/MatMul?
dense/BiasAdd/ReadVariableOpReadVariableOp%dense_biasadd_readvariableop_resource*
_output_shapes
:d*
dtype02
dense/BiasAdd/ReadVariableOp?
dense/BiasAddBiasAdddense/MatMul:product:0$dense/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????d2
dense/BiasAddg
	dense/EluEludense/BiasAdd:output:0*
T0*'
_output_shapes
:?????????d2
	dense/Elu?
dense_1/MatMul/ReadVariableOpReadVariableOp&dense_1_matmul_readvariableop_resource*
_output_shapes

:d
*
dtype02
dense_1/MatMul/ReadVariableOp?
dense_1/MatMulMatMuldense/Elu:activations:0%dense_1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????
2
dense_1/MatMul?
dense_1/BiasAdd/ReadVariableOpReadVariableOp'dense_1_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02 
dense_1/BiasAdd/ReadVariableOp?
dense_1/BiasAddBiasAdddense_1/MatMul:product:0&dense_1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????
2
dense_1/BiasAddm
dense_1/EluEludense_1/BiasAdd:output:0*
T0*'
_output_shapes
:?????????
2
dense_1/Elu?
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOp$dense_matmul_readvariableop_resource*
_output_shapes
:	?d*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOp?
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?d2!
dense/kernel/Regularizer/Square?
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const?
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/Sum?
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2 
dense/kernel/Regularizer/mul/x?
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/mul?
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_1_matmul_readvariableop_resource*
_output_shapes

:d
*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp?
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d
2#
!dense_1/kernel/Regularizer/Square?
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const?
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/Sum?
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2"
 dense_1/kernel/Regularizer/mul/x?
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/mulm
IdentityIdentitydense_1/Elu:activations:0*
T0*'
_output_shapes
:?????????
2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????:::::::::::::::::::::::::::::::\ X
4
_output_shapes"
 :??????????
 
_user_specified_nameinputs
?
?
D__inference_dense_2_layer_call_and_return_conditional_losses_1224532

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity??
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
MatMul?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2	
BiasAddd
IdentityIdentityBiasAdd:output:0*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*.
_input_shapes
:?????????:::O K
'
_output_shapes
:?????????
 
_user_specified_nameinputs
?
?
D__inference_dense_1_layer_call_and_return_conditional_losses_1223852

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity??
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:d
*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????
2
MatMul?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:
*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????
2	
BiasAddU
EluEluBiasAdd:output:0*
T0*'
_output_shapes
:?????????
2
Elu?
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:d
*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp?
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d
2#
!dense_1/kernel/Regularizer/Square?
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const?
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/Sum?
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2"
 dense_1/kernel/Regularizer/mul/x?
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/mule
IdentityIdentityElu:activations:0*
T0*'
_output_shapes
:?????????
2

Identity"
identityIdentity:output:0*.
_input_shapes
:?????????d:::O K
'
_output_shapes
:?????????d
 
_user_specified_nameinputs
?*
?
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_1227231

inputs
assignmovingavg_1227206
assignmovingavg_1_1227212)
%batchnorm_mul_readvariableop_resource%
!batchnorm_readvariableop_resource
identity??#AssignMovingAvg/AssignSubVariableOp?%AssignMovingAvg_1/AssignSubVariableOp?
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2 
moments/mean/reduction_indices?
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/mean?
moments/StopGradientStopGradientmoments/mean:output:0*
T0**
_output_shapes
:2
moments/StopGradient?
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*3
_output_shapes!
:?????????2
moments/SquaredDifference?
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2$
"moments/variance/reduction_indices?
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/variance?
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze?
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze_1?
AssignMovingAvg/decayConst**
_class 
loc:@AssignMovingAvg/1227206*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg/decay?
AssignMovingAvg/ReadVariableOpReadVariableOpassignmovingavg_1227206*
_output_shapes
:*
dtype02 
AssignMovingAvg/ReadVariableOp?
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0**
_class 
loc:@AssignMovingAvg/1227206*
_output_shapes
:2
AssignMovingAvg/sub?
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0**
_class 
loc:@AssignMovingAvg/1227206*
_output_shapes
:2
AssignMovingAvg/mul?
#AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1227206AssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp**
_class 
loc:@AssignMovingAvg/1227206*
_output_shapes
 *
dtype02%
#AssignMovingAvg/AssignSubVariableOp?
AssignMovingAvg_1/decayConst*,
_class"
 loc:@AssignMovingAvg_1/1227212*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg_1/decay?
 AssignMovingAvg_1/ReadVariableOpReadVariableOpassignmovingavg_1_1227212*
_output_shapes
:*
dtype02"
 AssignMovingAvg_1/ReadVariableOp?
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1227212*
_output_shapes
:2
AssignMovingAvg_1/sub?
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1227212*
_output_shapes
:2
AssignMovingAvg_1/mul?
%AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1_1227212AssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*,
_class"
 loc:@AssignMovingAvg_1/1227212*
_output_shapes
 *
dtype02'
%AssignMovingAvg_1/AssignSubVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/mul_1{
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp?
batchnorm/subSub batchnorm/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0$^AssignMovingAvg/AssignSubVariableOp&^AssignMovingAvg_1/AssignSubVariableOp*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????::::2J
#AssignMovingAvg/AssignSubVariableOp#AssignMovingAvg/AssignSubVariableOp2N
%AssignMovingAvg_1/AssignSubVariableOp%AssignMovingAvg_1/AssignSubVariableOp:[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?
}
(__inference_conv3d_layer_call_fn_1226929

inputs
unknown
	unknown_0
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *L
fGRE
C__inference_conv3d_layer_call_and_return_conditional_losses_12227782
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*;
_input_shapes*
(:??????????::22
StatefulPartitionedCallStatefulPartitionedCall:\ X
4
_output_shapes"
 :??????????
 
_user_specified_nameinputs
?
m
C__inference_lambda_layer_call_and_return_conditional_losses_1224484

inputs
inputs_1
identityU
subSubinputsinputs_1*
T0*'
_output_shapes
:?????????
2
subL
AbsAbssub:z:0*
T0*'
_output_shapes
:?????????
2
Abs[
IdentityIdentityAbs:y:0*
T0*'
_output_shapes
:?????????
2

Identity"
identityIdentity:output:0*9
_input_shapes(
&:?????????
:?????????
:O K
'
_output_shapes
:?????????

 
_user_specified_nameinputs:OK
'
_output_shapes
:?????????

 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_1227537

inputs%
!batchnorm_readvariableop_resource)
%batchnorm_mul_readvariableop_resource'
#batchnorm_readvariableop_1_resource'
#batchnorm_readvariableop_2_resource
identity??
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/mul_1?
batchnorm/ReadVariableOp_1ReadVariableOp#batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_1?
batchnorm/mul_2Mul"batchnorm/ReadVariableOp_1:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOp_2ReadVariableOp#batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_2?
batchnorm/subSub"batchnorm/ReadVariableOp_2:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????:::::v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
?
b
D__inference_dropout_layer_call_and_return_conditional_losses_1223271

inputs

identity_1[
IdentityIdentityinputs*
T0*(
_output_shapes
:??????????
2

Identityj

Identity_1IdentityIdentity:output:0*
T0*(
_output_shapes
:??????????
2

Identity_1"!

identity_1Identity_1:output:0*'
_input_shapes
:??????????
:P L
(
_output_shapes
:??????????

 
_user_specified_nameinputs
?
?
,__inference_sequential_layer_call_fn_1226369

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16

unknown_17

unknown_18

unknown_19

unknown_20

unknown_21

unknown_22

unknown_23

unknown_24

unknown_25

unknown_26

unknown_27

unknown_28
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19
unknown_20
unknown_21
unknown_22
unknown_23
unknown_24
unknown_25
unknown_26
unknown_27
unknown_28**
Tin#
!2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
*@
_read_only_resource_inputs"
 	
*0
config_proto 

CPU

GPU2*0J 8? *P
fKRI
G__inference_sequential_layer_call_and_return_conditional_losses_12241862
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????
2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????::::::::::::::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:\ X
4
_output_shapes"
 :??????????
 
_user_specified_nameinputs
?
O
3__inference_average_pooling3d_layer_call_fn_1222764

inputs
identity?
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *W
_output_shapesE
C:A?????????????????????????????????????????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *W
fRRP
N__inference_average_pooling3d_layer_call_and_return_conditional_losses_12227582
PartitionedCall?
IdentityIdentityPartitionedCall:output:0*
T0*W
_output_shapesE
C:A?????????????????????????????????????????????2

Identity"
identityIdentity:output:0*V
_input_shapesE
C:A?????????????????????????????????????????????: {
W
_output_shapesE
C:A?????????????????????????????????????????????
 
_user_specified_nameinputs
?
c
D__inference_dropout_layer_call_and_return_conditional_losses_1227668

inputs
identity?c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *UU??2
dropout/Constt
dropout/MulMulinputsdropout/Const:output:0*
T0*(
_output_shapes
:??????????
2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape?
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*(
_output_shapes
:??????????
*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2
dropout/GreaterEqual/y?
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*(
_output_shapes
:??????????
2
dropout/GreaterEqual?
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*(
_output_shapes
:??????????
2
dropout/Cast{
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*(
_output_shapes
:??????????
2
dropout/Mul_1f
IdentityIdentitydropout/Mul_1:z:0*
T0*(
_output_shapes
:??????????
2

Identity"
identityIdentity:output:0*'
_input_shapes
:??????????
:P L
(
_output_shapes
:??????????

 
_user_specified_nameinputs
?
?
%__inference_signature_wrapper_1225148
input_1
input_2
input_3
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16

unknown_17

unknown_18

unknown_19

unknown_20

unknown_21

unknown_22

unknown_23

unknown_24

unknown_25

unknown_26

unknown_27

unknown_28

unknown_29

unknown_30
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinput_1input_2input_3unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19
unknown_20
unknown_21
unknown_22
unknown_23
unknown_24
unknown_25
unknown_26
unknown_27
unknown_28
unknown_29
unknown_30*.
Tin'
%2#*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????*B
_read_only_resource_inputs$
" 	
 !"*0
config_proto 

CPU

GPU2*0J 8? *+
f&R$
"__inference__wrapped_model_12221922
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????:??????????:?????????
::::::::::::::::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:] Y
4
_output_shapes"
 :??????????
!
_user_specified_name	input_1:]Y
4
_output_shapes"
 :??????????
!
_user_specified_name	input_2:PL
'
_output_shapes
:?????????

!
_user_specified_name	input_3
?
?
.__inference_functional_3_layer_call_fn_1225821
inputs_0
inputs_1
inputs_2
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16

unknown_17

unknown_18

unknown_19

unknown_20

unknown_21

unknown_22

unknown_23

unknown_24

unknown_25

unknown_26

unknown_27

unknown_28

unknown_29

unknown_30
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputs_0inputs_1inputs_2unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19
unknown_20
unknown_21
unknown_22
unknown_23
unknown_24
unknown_25
unknown_26
unknown_27
unknown_28
unknown_29
unknown_30*.
Tin'
%2#*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????*:
_read_only_resource_inputs
	
 !"*0
config_proto 

CPU

GPU2*0J 8? *R
fMRK
I__inference_functional_3_layer_call_and_return_conditional_losses_12248002
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????:??????????:?????????
::::::::::::::::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:^ Z
4
_output_shapes"
 :??????????
"
_user_specified_name
inputs/0:^Z
4
_output_shapes"
 :??????????
"
_user_specified_name
inputs/1:QM
'
_output_shapes
:?????????

"
_user_specified_name
inputs/2
?,
?
G__inference_sequential_layer_call_and_return_conditional_losses_1223960
functional_1_input
functional_1_1223884
functional_1_1223886
functional_1_1223888
functional_1_1223890
functional_1_1223892
functional_1_1223894
functional_1_1223896
functional_1_1223898
functional_1_1223900
functional_1_1223902
functional_1_1223904
functional_1_1223906
functional_1_1223908
functional_1_1223910
functional_1_1223912
functional_1_1223914
functional_1_1223916
functional_1_1223918
functional_1_1223920
functional_1_1223922
functional_1_1223924
functional_1_1223926
functional_1_1223928
functional_1_1223930
functional_1_1223932
functional_1_1223934
dense_1223937
dense_1223939
dense_1_1223942
dense_1_1223944
identity??dense/StatefulPartitionedCall?dense_1/StatefulPartitionedCall?$functional_1/StatefulPartitionedCall?
$functional_1/StatefulPartitionedCallStatefulPartitionedCallfunctional_1_inputfunctional_1_1223884functional_1_1223886functional_1_1223888functional_1_1223890functional_1_1223892functional_1_1223894functional_1_1223896functional_1_1223898functional_1_1223900functional_1_1223902functional_1_1223904functional_1_1223906functional_1_1223908functional_1_1223910functional_1_1223912functional_1_1223914functional_1_1223916functional_1_1223918functional_1_1223920functional_1_1223922functional_1_1223924functional_1_1223926functional_1_1223928functional_1_1223930functional_1_1223932functional_1_1223934*&
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*<
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8? *R
fMRK
I__inference_functional_1_layer_call_and_return_conditional_losses_12235762&
$functional_1/StatefulPartitionedCall?
dense/StatefulPartitionedCallStatefulPartitionedCall-functional_1/StatefulPartitionedCall:output:0dense_1223937dense_1223939*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????d*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *K
fFRD
B__inference_dense_layer_call_and_return_conditional_losses_12238192
dense/StatefulPartitionedCall?
dense_1/StatefulPartitionedCallStatefulPartitionedCall&dense/StatefulPartitionedCall:output:0dense_1_1223942dense_1_1223944*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_1_layer_call_and_return_conditional_losses_12238522!
dense_1/StatefulPartitionedCall?
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_1223937*
_output_shapes
:	?d*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOp?
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?d2!
dense/kernel/Regularizer/Square?
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const?
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/Sum?
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2 
dense/kernel/Regularizer/mul/x?
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/mul?
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_1223942*
_output_shapes

:d
*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp?
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d
2#
!dense_1/kernel/Regularizer/Square?
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const?
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/Sum?
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2"
 dense_1/kernel/Regularizer/mul/x?
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/mul?
IdentityIdentity(dense_1/StatefulPartitionedCall:output:0^dense/StatefulPartitionedCall ^dense_1/StatefulPartitionedCall%^functional_1/StatefulPartitionedCall*
T0*'
_output_shapes
:?????????
2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????::::::::::::::::::::::::::::::2>
dense/StatefulPartitionedCalldense/StatefulPartitionedCall2B
dense_1/StatefulPartitionedCalldense_1/StatefulPartitionedCall2L
$functional_1/StatefulPartitionedCall$functional_1/StatefulPartitionedCall:h d
4
_output_shapes"
 :??????????
,
_user_specified_namefunctional_1_input
?
?
.__inference_functional_3_layer_call_fn_1225892
inputs_0
inputs_1
inputs_2
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16

unknown_17

unknown_18

unknown_19

unknown_20

unknown_21

unknown_22

unknown_23

unknown_24

unknown_25

unknown_26

unknown_27

unknown_28

unknown_29

unknown_30
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputs_0inputs_1inputs_2unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19
unknown_20
unknown_21
unknown_22
unknown_23
unknown_24
unknown_25
unknown_26
unknown_27
unknown_28
unknown_29
unknown_30*.
Tin'
%2#*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????*B
_read_only_resource_inputs$
" 	
 !"*0
config_proto 

CPU

GPU2*0J 8? *R
fMRK
I__inference_functional_3_layer_call_and_return_conditional_losses_12249882
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????:??????????:?????????
::::::::::::::::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:^ Z
4
_output_shapes"
 :??????????
"
_user_specified_name
inputs/0:^Z
4
_output_shapes"
 :??????????
"
_user_specified_name
inputs/1:QM
'
_output_shapes
:?????????

"
_user_specified_name
inputs/2
?
?
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_1223085

inputs%
!batchnorm_readvariableop_resource)
%batchnorm_mul_readvariableop_resource'
#batchnorm_readvariableop_1_resource'
#batchnorm_readvariableop_2_resource
identity??
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/mul_1?
batchnorm/ReadVariableOp_1ReadVariableOp#batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_1?
batchnorm/mul_2Mul"batchnorm/ReadVariableOp_1:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOp_2ReadVariableOp#batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_2?
batchnorm/subSub"batchnorm/ReadVariableOp_2:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/add_1s
IdentityIdentitybatchnorm/add_1:z:0*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????:::::[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?
~
)__inference_dense_2_layer_call_fn_1226427

inputs
unknown
	unknown_0
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_2_layer_call_and_return_conditional_losses_12245322
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*.
_input_shapes
:?????????::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:?????????
 
_user_specified_nameinputs
?
c
D__inference_dropout_layer_call_and_return_conditional_losses_1223266

inputs
identity?c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *UU??2
dropout/Constt
dropout/MulMulinputsdropout/Const:output:0*
T0*(
_output_shapes
:??????????
2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape?
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*(
_output_shapes
:??????????
*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2
dropout/GreaterEqual/y?
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*(
_output_shapes
:??????????
2
dropout/GreaterEqual?
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*(
_output_shapes
:??????????
2
dropout/Cast{
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*(
_output_shapes
:??????????
2
dropout/Mul_1f
IdentityIdentitydropout/Mul_1:z:0*
T0*(
_output_shapes
:??????????
2

Identity"
identityIdentity:output:0*'
_input_shapes
:??????????
:P L
(
_output_shapes
:??????????

 
_user_specified_nameinputs
?
?
C__inference_layer1_layer_call_and_return_conditional_losses_1223295

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity??
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource* 
_output_shapes
:
?
?*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2
MatMul?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:?*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2	
BiasAddV
EluEluBiasAdd:output:0*
T0*(
_output_shapes
:??????????2
Eluf
IdentityIdentityElu:activations:0*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*/
_input_shapes
:??????????
:::P L
(
_output_shapes
:??????????

 
_user_specified_nameinputs
?
?
7__inference_batch_normalization_1_layer_call_fn_1227182

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *N
_output_shapes<
::8????????????????????????????????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_12224282
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????::::22
StatefulPartitionedCallStatefulPartitionedCall:v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
?+
?
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_1227149

inputs
assignmovingavg_1227124
assignmovingavg_1_1227130)
%batchnorm_mul_readvariableop_resource%
!batchnorm_readvariableop_resource
identity??#AssignMovingAvg/AssignSubVariableOp?%AssignMovingAvg_1/AssignSubVariableOp?
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2 
moments/mean/reduction_indices?
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/mean?
moments/StopGradientStopGradientmoments/mean:output:0*
T0**
_output_shapes
:2
moments/StopGradient?
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
moments/SquaredDifference?
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2$
"moments/variance/reduction_indices?
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/variance?
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze?
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze_1?
AssignMovingAvg/decayConst**
_class 
loc:@AssignMovingAvg/1227124*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg/decay?
AssignMovingAvg/ReadVariableOpReadVariableOpassignmovingavg_1227124*
_output_shapes
:*
dtype02 
AssignMovingAvg/ReadVariableOp?
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0**
_class 
loc:@AssignMovingAvg/1227124*
_output_shapes
:2
AssignMovingAvg/sub?
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0**
_class 
loc:@AssignMovingAvg/1227124*
_output_shapes
:2
AssignMovingAvg/mul?
#AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1227124AssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp**
_class 
loc:@AssignMovingAvg/1227124*
_output_shapes
 *
dtype02%
#AssignMovingAvg/AssignSubVariableOp?
AssignMovingAvg_1/decayConst*,
_class"
 loc:@AssignMovingAvg_1/1227130*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg_1/decay?
 AssignMovingAvg_1/ReadVariableOpReadVariableOpassignmovingavg_1_1227130*
_output_shapes
:*
dtype02"
 AssignMovingAvg_1/ReadVariableOp?
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1227130*
_output_shapes
:2
AssignMovingAvg_1/sub?
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1227130*
_output_shapes
:2
AssignMovingAvg_1/mul?
%AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1_1227130AssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*,
_class"
 loc:@AssignMovingAvg_1/1227130*
_output_shapes
 *
dtype02'
%AssignMovingAvg_1/AssignSubVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/mul_1{
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp?
batchnorm/subSub batchnorm/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0$^AssignMovingAvg/AssignSubVariableOp&^AssignMovingAvg_1/AssignSubVariableOp*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????::::2J
#AssignMovingAvg/AssignSubVariableOp#AssignMovingAvg/AssignSubVariableOp2N
%AssignMovingAvg_1/AssignSubVariableOp%AssignMovingAvg_1/AssignSubVariableOp:v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
??
?0
#__inference__traced_restore_1228260
file_prefix#
assignvariableop_dense_2_kernel#
assignvariableop_1_dense_2_bias 
assignvariableop_2_adam_iter"
assignvariableop_3_adam_beta_1"
assignvariableop_4_adam_beta_2!
assignvariableop_5_adam_decay)
%assignvariableop_6_adam_learning_rate$
 assignvariableop_7_conv3d_kernel"
assignvariableop_8_conv3d_bias0
,assignvariableop_9_batch_normalization_gamma0
,assignvariableop_10_batch_normalization_beta'
#assignvariableop_11_conv3d_1_kernel%
!assignvariableop_12_conv3d_1_bias3
/assignvariableop_13_batch_normalization_1_gamma2
.assignvariableop_14_batch_normalization_1_beta'
#assignvariableop_15_conv3d_2_kernel%
!assignvariableop_16_conv3d_2_bias3
/assignvariableop_17_batch_normalization_2_gamma2
.assignvariableop_18_batch_normalization_2_beta'
#assignvariableop_19_conv3d_3_kernel%
!assignvariableop_20_conv3d_3_bias3
/assignvariableop_21_batch_normalization_3_gamma2
.assignvariableop_22_batch_normalization_3_beta%
!assignvariableop_23_layer1_kernel#
assignvariableop_24_layer1_bias$
 assignvariableop_25_dense_kernel"
assignvariableop_26_dense_bias&
"assignvariableop_27_dense_1_kernel$
 assignvariableop_28_dense_1_bias7
3assignvariableop_29_batch_normalization_moving_mean;
7assignvariableop_30_batch_normalization_moving_variance9
5assignvariableop_31_batch_normalization_1_moving_mean=
9assignvariableop_32_batch_normalization_1_moving_variance9
5assignvariableop_33_batch_normalization_2_moving_mean=
9assignvariableop_34_batch_normalization_2_moving_variance9
5assignvariableop_35_batch_normalization_3_moving_mean=
9assignvariableop_36_batch_normalization_3_moving_variance
assignvariableop_37_total
assignvariableop_38_count-
)assignvariableop_39_adam_dense_2_kernel_m+
'assignvariableop_40_adam_dense_2_bias_m,
(assignvariableop_41_adam_conv3d_kernel_m*
&assignvariableop_42_adam_conv3d_bias_m8
4assignvariableop_43_adam_batch_normalization_gamma_m7
3assignvariableop_44_adam_batch_normalization_beta_m.
*assignvariableop_45_adam_conv3d_1_kernel_m,
(assignvariableop_46_adam_conv3d_1_bias_m:
6assignvariableop_47_adam_batch_normalization_1_gamma_m9
5assignvariableop_48_adam_batch_normalization_1_beta_m.
*assignvariableop_49_adam_conv3d_2_kernel_m,
(assignvariableop_50_adam_conv3d_2_bias_m:
6assignvariableop_51_adam_batch_normalization_2_gamma_m9
5assignvariableop_52_adam_batch_normalization_2_beta_m.
*assignvariableop_53_adam_conv3d_3_kernel_m,
(assignvariableop_54_adam_conv3d_3_bias_m:
6assignvariableop_55_adam_batch_normalization_3_gamma_m9
5assignvariableop_56_adam_batch_normalization_3_beta_m,
(assignvariableop_57_adam_layer1_kernel_m*
&assignvariableop_58_adam_layer1_bias_m+
'assignvariableop_59_adam_dense_kernel_m)
%assignvariableop_60_adam_dense_bias_m-
)assignvariableop_61_adam_dense_1_kernel_m+
'assignvariableop_62_adam_dense_1_bias_m-
)assignvariableop_63_adam_dense_2_kernel_v+
'assignvariableop_64_adam_dense_2_bias_v,
(assignvariableop_65_adam_conv3d_kernel_v*
&assignvariableop_66_adam_conv3d_bias_v8
4assignvariableop_67_adam_batch_normalization_gamma_v7
3assignvariableop_68_adam_batch_normalization_beta_v.
*assignvariableop_69_adam_conv3d_1_kernel_v,
(assignvariableop_70_adam_conv3d_1_bias_v:
6assignvariableop_71_adam_batch_normalization_1_gamma_v9
5assignvariableop_72_adam_batch_normalization_1_beta_v.
*assignvariableop_73_adam_conv3d_2_kernel_v,
(assignvariableop_74_adam_conv3d_2_bias_v:
6assignvariableop_75_adam_batch_normalization_2_gamma_v9
5assignvariableop_76_adam_batch_normalization_2_beta_v.
*assignvariableop_77_adam_conv3d_3_kernel_v,
(assignvariableop_78_adam_conv3d_3_bias_v:
6assignvariableop_79_adam_batch_normalization_3_gamma_v9
5assignvariableop_80_adam_batch_normalization_3_beta_v,
(assignvariableop_81_adam_layer1_kernel_v*
&assignvariableop_82_adam_layer1_bias_v+
'assignvariableop_83_adam_dense_kernel_v)
%assignvariableop_84_adam_dense_bias_v-
)assignvariableop_85_adam_dense_1_kernel_v+
'assignvariableop_86_adam_dense_1_bias_v
identity_88??AssignVariableOp?AssignVariableOp_1?AssignVariableOp_10?AssignVariableOp_11?AssignVariableOp_12?AssignVariableOp_13?AssignVariableOp_14?AssignVariableOp_15?AssignVariableOp_16?AssignVariableOp_17?AssignVariableOp_18?AssignVariableOp_19?AssignVariableOp_2?AssignVariableOp_20?AssignVariableOp_21?AssignVariableOp_22?AssignVariableOp_23?AssignVariableOp_24?AssignVariableOp_25?AssignVariableOp_26?AssignVariableOp_27?AssignVariableOp_28?AssignVariableOp_29?AssignVariableOp_3?AssignVariableOp_30?AssignVariableOp_31?AssignVariableOp_32?AssignVariableOp_33?AssignVariableOp_34?AssignVariableOp_35?AssignVariableOp_36?AssignVariableOp_37?AssignVariableOp_38?AssignVariableOp_39?AssignVariableOp_4?AssignVariableOp_40?AssignVariableOp_41?AssignVariableOp_42?AssignVariableOp_43?AssignVariableOp_44?AssignVariableOp_45?AssignVariableOp_46?AssignVariableOp_47?AssignVariableOp_48?AssignVariableOp_49?AssignVariableOp_5?AssignVariableOp_50?AssignVariableOp_51?AssignVariableOp_52?AssignVariableOp_53?AssignVariableOp_54?AssignVariableOp_55?AssignVariableOp_56?AssignVariableOp_57?AssignVariableOp_58?AssignVariableOp_59?AssignVariableOp_6?AssignVariableOp_60?AssignVariableOp_61?AssignVariableOp_62?AssignVariableOp_63?AssignVariableOp_64?AssignVariableOp_65?AssignVariableOp_66?AssignVariableOp_67?AssignVariableOp_68?AssignVariableOp_69?AssignVariableOp_7?AssignVariableOp_70?AssignVariableOp_71?AssignVariableOp_72?AssignVariableOp_73?AssignVariableOp_74?AssignVariableOp_75?AssignVariableOp_76?AssignVariableOp_77?AssignVariableOp_78?AssignVariableOp_79?AssignVariableOp_8?AssignVariableOp_80?AssignVariableOp_81?AssignVariableOp_82?AssignVariableOp_83?AssignVariableOp_84?AssignVariableOp_85?AssignVariableOp_86?AssignVariableOp_9?-
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:X*
dtype0*?,
value?,B?,XB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/0/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/1/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/2/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/3/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/4/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/5/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/6/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/7/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/8/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/9/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/10/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/11/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/12/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/13/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/14/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/15/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/16/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/17/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/18/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/19/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/20/.ATTRIBUTES/VARIABLE_VALUEB1trainable_variables/21/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB'variables/10/.ATTRIBUTES/VARIABLE_VALUEB'variables/11/.ATTRIBUTES/VARIABLE_VALUEB'variables/16/.ATTRIBUTES/VARIABLE_VALUEB'variables/17/.ATTRIBUTES/VARIABLE_VALUEB'variables/22/.ATTRIBUTES/VARIABLE_VALUEB'variables/23/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/0/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/1/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/4/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/5/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/7/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/8/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/9/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/10/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/11/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/12/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/13/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/14/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/15/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/16/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/17/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/18/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/19/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/20/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/21/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-1/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-1/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/0/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/1/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/4/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/5/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/7/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/8/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/9/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/10/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/11/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/12/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/13/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/14/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/15/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/16/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/17/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/18/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/19/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/20/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBMtrainable_variables/21/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH2
RestoreV2/tensor_names?
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:X*
dtype0*?
value?B?XB B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 2
RestoreV2/shape_and_slices?
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*?
_output_shapes?
?::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*f
dtypes\
Z2X	2
	RestoreV2g
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0*
_output_shapes
:2

Identity?
AssignVariableOpAssignVariableOpassignvariableop_dense_2_kernelIdentity:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOpk

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:2

Identity_1?
AssignVariableOp_1AssignVariableOpassignvariableop_1_dense_2_biasIdentity_1:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_1k

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0	*
_output_shapes
:2

Identity_2?
AssignVariableOp_2AssignVariableOpassignvariableop_2_adam_iterIdentity_2:output:0"/device:CPU:0*
_output_shapes
 *
dtype0	2
AssignVariableOp_2k

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:2

Identity_3?
AssignVariableOp_3AssignVariableOpassignvariableop_3_adam_beta_1Identity_3:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_3k

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:2

Identity_4?
AssignVariableOp_4AssignVariableOpassignvariableop_4_adam_beta_2Identity_4:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_4k

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:2

Identity_5?
AssignVariableOp_5AssignVariableOpassignvariableop_5_adam_decayIdentity_5:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_5k

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:2

Identity_6?
AssignVariableOp_6AssignVariableOp%assignvariableop_6_adam_learning_rateIdentity_6:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_6k

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:2

Identity_7?
AssignVariableOp_7AssignVariableOp assignvariableop_7_conv3d_kernelIdentity_7:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_7k

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:2

Identity_8?
AssignVariableOp_8AssignVariableOpassignvariableop_8_conv3d_biasIdentity_8:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_8k

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:2

Identity_9?
AssignVariableOp_9AssignVariableOp,assignvariableop_9_batch_normalization_gammaIdentity_9:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_9n
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:2
Identity_10?
AssignVariableOp_10AssignVariableOp,assignvariableop_10_batch_normalization_betaIdentity_10:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_10n
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:2
Identity_11?
AssignVariableOp_11AssignVariableOp#assignvariableop_11_conv3d_1_kernelIdentity_11:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_11n
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0*
_output_shapes
:2
Identity_12?
AssignVariableOp_12AssignVariableOp!assignvariableop_12_conv3d_1_biasIdentity_12:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_12n
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:2
Identity_13?
AssignVariableOp_13AssignVariableOp/assignvariableop_13_batch_normalization_1_gammaIdentity_13:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_13n
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0*
_output_shapes
:2
Identity_14?
AssignVariableOp_14AssignVariableOp.assignvariableop_14_batch_normalization_1_betaIdentity_14:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_14n
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:2
Identity_15?
AssignVariableOp_15AssignVariableOp#assignvariableop_15_conv3d_2_kernelIdentity_15:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_15n
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:2
Identity_16?
AssignVariableOp_16AssignVariableOp!assignvariableop_16_conv3d_2_biasIdentity_16:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_16n
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:2
Identity_17?
AssignVariableOp_17AssignVariableOp/assignvariableop_17_batch_normalization_2_gammaIdentity_17:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_17n
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:2
Identity_18?
AssignVariableOp_18AssignVariableOp.assignvariableop_18_batch_normalization_2_betaIdentity_18:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_18n
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:2
Identity_19?
AssignVariableOp_19AssignVariableOp#assignvariableop_19_conv3d_3_kernelIdentity_19:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_19n
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0*
_output_shapes
:2
Identity_20?
AssignVariableOp_20AssignVariableOp!assignvariableop_20_conv3d_3_biasIdentity_20:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_20n
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:2
Identity_21?
AssignVariableOp_21AssignVariableOp/assignvariableop_21_batch_normalization_3_gammaIdentity_21:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_21n
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:2
Identity_22?
AssignVariableOp_22AssignVariableOp.assignvariableop_22_batch_normalization_3_betaIdentity_22:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_22n
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:2
Identity_23?
AssignVariableOp_23AssignVariableOp!assignvariableop_23_layer1_kernelIdentity_23:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_23n
Identity_24IdentityRestoreV2:tensors:24"/device:CPU:0*
T0*
_output_shapes
:2
Identity_24?
AssignVariableOp_24AssignVariableOpassignvariableop_24_layer1_biasIdentity_24:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_24n
Identity_25IdentityRestoreV2:tensors:25"/device:CPU:0*
T0*
_output_shapes
:2
Identity_25?
AssignVariableOp_25AssignVariableOp assignvariableop_25_dense_kernelIdentity_25:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_25n
Identity_26IdentityRestoreV2:tensors:26"/device:CPU:0*
T0*
_output_shapes
:2
Identity_26?
AssignVariableOp_26AssignVariableOpassignvariableop_26_dense_biasIdentity_26:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_26n
Identity_27IdentityRestoreV2:tensors:27"/device:CPU:0*
T0*
_output_shapes
:2
Identity_27?
AssignVariableOp_27AssignVariableOp"assignvariableop_27_dense_1_kernelIdentity_27:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_27n
Identity_28IdentityRestoreV2:tensors:28"/device:CPU:0*
T0*
_output_shapes
:2
Identity_28?
AssignVariableOp_28AssignVariableOp assignvariableop_28_dense_1_biasIdentity_28:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_28n
Identity_29IdentityRestoreV2:tensors:29"/device:CPU:0*
T0*
_output_shapes
:2
Identity_29?
AssignVariableOp_29AssignVariableOp3assignvariableop_29_batch_normalization_moving_meanIdentity_29:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_29n
Identity_30IdentityRestoreV2:tensors:30"/device:CPU:0*
T0*
_output_shapes
:2
Identity_30?
AssignVariableOp_30AssignVariableOp7assignvariableop_30_batch_normalization_moving_varianceIdentity_30:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_30n
Identity_31IdentityRestoreV2:tensors:31"/device:CPU:0*
T0*
_output_shapes
:2
Identity_31?
AssignVariableOp_31AssignVariableOp5assignvariableop_31_batch_normalization_1_moving_meanIdentity_31:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_31n
Identity_32IdentityRestoreV2:tensors:32"/device:CPU:0*
T0*
_output_shapes
:2
Identity_32?
AssignVariableOp_32AssignVariableOp9assignvariableop_32_batch_normalization_1_moving_varianceIdentity_32:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_32n
Identity_33IdentityRestoreV2:tensors:33"/device:CPU:0*
T0*
_output_shapes
:2
Identity_33?
AssignVariableOp_33AssignVariableOp5assignvariableop_33_batch_normalization_2_moving_meanIdentity_33:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_33n
Identity_34IdentityRestoreV2:tensors:34"/device:CPU:0*
T0*
_output_shapes
:2
Identity_34?
AssignVariableOp_34AssignVariableOp9assignvariableop_34_batch_normalization_2_moving_varianceIdentity_34:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_34n
Identity_35IdentityRestoreV2:tensors:35"/device:CPU:0*
T0*
_output_shapes
:2
Identity_35?
AssignVariableOp_35AssignVariableOp5assignvariableop_35_batch_normalization_3_moving_meanIdentity_35:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_35n
Identity_36IdentityRestoreV2:tensors:36"/device:CPU:0*
T0*
_output_shapes
:2
Identity_36?
AssignVariableOp_36AssignVariableOp9assignvariableop_36_batch_normalization_3_moving_varianceIdentity_36:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_36n
Identity_37IdentityRestoreV2:tensors:37"/device:CPU:0*
T0*
_output_shapes
:2
Identity_37?
AssignVariableOp_37AssignVariableOpassignvariableop_37_totalIdentity_37:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_37n
Identity_38IdentityRestoreV2:tensors:38"/device:CPU:0*
T0*
_output_shapes
:2
Identity_38?
AssignVariableOp_38AssignVariableOpassignvariableop_38_countIdentity_38:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_38n
Identity_39IdentityRestoreV2:tensors:39"/device:CPU:0*
T0*
_output_shapes
:2
Identity_39?
AssignVariableOp_39AssignVariableOp)assignvariableop_39_adam_dense_2_kernel_mIdentity_39:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_39n
Identity_40IdentityRestoreV2:tensors:40"/device:CPU:0*
T0*
_output_shapes
:2
Identity_40?
AssignVariableOp_40AssignVariableOp'assignvariableop_40_adam_dense_2_bias_mIdentity_40:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_40n
Identity_41IdentityRestoreV2:tensors:41"/device:CPU:0*
T0*
_output_shapes
:2
Identity_41?
AssignVariableOp_41AssignVariableOp(assignvariableop_41_adam_conv3d_kernel_mIdentity_41:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_41n
Identity_42IdentityRestoreV2:tensors:42"/device:CPU:0*
T0*
_output_shapes
:2
Identity_42?
AssignVariableOp_42AssignVariableOp&assignvariableop_42_adam_conv3d_bias_mIdentity_42:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_42n
Identity_43IdentityRestoreV2:tensors:43"/device:CPU:0*
T0*
_output_shapes
:2
Identity_43?
AssignVariableOp_43AssignVariableOp4assignvariableop_43_adam_batch_normalization_gamma_mIdentity_43:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_43n
Identity_44IdentityRestoreV2:tensors:44"/device:CPU:0*
T0*
_output_shapes
:2
Identity_44?
AssignVariableOp_44AssignVariableOp3assignvariableop_44_adam_batch_normalization_beta_mIdentity_44:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_44n
Identity_45IdentityRestoreV2:tensors:45"/device:CPU:0*
T0*
_output_shapes
:2
Identity_45?
AssignVariableOp_45AssignVariableOp*assignvariableop_45_adam_conv3d_1_kernel_mIdentity_45:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_45n
Identity_46IdentityRestoreV2:tensors:46"/device:CPU:0*
T0*
_output_shapes
:2
Identity_46?
AssignVariableOp_46AssignVariableOp(assignvariableop_46_adam_conv3d_1_bias_mIdentity_46:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_46n
Identity_47IdentityRestoreV2:tensors:47"/device:CPU:0*
T0*
_output_shapes
:2
Identity_47?
AssignVariableOp_47AssignVariableOp6assignvariableop_47_adam_batch_normalization_1_gamma_mIdentity_47:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_47n
Identity_48IdentityRestoreV2:tensors:48"/device:CPU:0*
T0*
_output_shapes
:2
Identity_48?
AssignVariableOp_48AssignVariableOp5assignvariableop_48_adam_batch_normalization_1_beta_mIdentity_48:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_48n
Identity_49IdentityRestoreV2:tensors:49"/device:CPU:0*
T0*
_output_shapes
:2
Identity_49?
AssignVariableOp_49AssignVariableOp*assignvariableop_49_adam_conv3d_2_kernel_mIdentity_49:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_49n
Identity_50IdentityRestoreV2:tensors:50"/device:CPU:0*
T0*
_output_shapes
:2
Identity_50?
AssignVariableOp_50AssignVariableOp(assignvariableop_50_adam_conv3d_2_bias_mIdentity_50:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_50n
Identity_51IdentityRestoreV2:tensors:51"/device:CPU:0*
T0*
_output_shapes
:2
Identity_51?
AssignVariableOp_51AssignVariableOp6assignvariableop_51_adam_batch_normalization_2_gamma_mIdentity_51:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_51n
Identity_52IdentityRestoreV2:tensors:52"/device:CPU:0*
T0*
_output_shapes
:2
Identity_52?
AssignVariableOp_52AssignVariableOp5assignvariableop_52_adam_batch_normalization_2_beta_mIdentity_52:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_52n
Identity_53IdentityRestoreV2:tensors:53"/device:CPU:0*
T0*
_output_shapes
:2
Identity_53?
AssignVariableOp_53AssignVariableOp*assignvariableop_53_adam_conv3d_3_kernel_mIdentity_53:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_53n
Identity_54IdentityRestoreV2:tensors:54"/device:CPU:0*
T0*
_output_shapes
:2
Identity_54?
AssignVariableOp_54AssignVariableOp(assignvariableop_54_adam_conv3d_3_bias_mIdentity_54:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_54n
Identity_55IdentityRestoreV2:tensors:55"/device:CPU:0*
T0*
_output_shapes
:2
Identity_55?
AssignVariableOp_55AssignVariableOp6assignvariableop_55_adam_batch_normalization_3_gamma_mIdentity_55:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_55n
Identity_56IdentityRestoreV2:tensors:56"/device:CPU:0*
T0*
_output_shapes
:2
Identity_56?
AssignVariableOp_56AssignVariableOp5assignvariableop_56_adam_batch_normalization_3_beta_mIdentity_56:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_56n
Identity_57IdentityRestoreV2:tensors:57"/device:CPU:0*
T0*
_output_shapes
:2
Identity_57?
AssignVariableOp_57AssignVariableOp(assignvariableop_57_adam_layer1_kernel_mIdentity_57:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_57n
Identity_58IdentityRestoreV2:tensors:58"/device:CPU:0*
T0*
_output_shapes
:2
Identity_58?
AssignVariableOp_58AssignVariableOp&assignvariableop_58_adam_layer1_bias_mIdentity_58:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_58n
Identity_59IdentityRestoreV2:tensors:59"/device:CPU:0*
T0*
_output_shapes
:2
Identity_59?
AssignVariableOp_59AssignVariableOp'assignvariableop_59_adam_dense_kernel_mIdentity_59:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_59n
Identity_60IdentityRestoreV2:tensors:60"/device:CPU:0*
T0*
_output_shapes
:2
Identity_60?
AssignVariableOp_60AssignVariableOp%assignvariableop_60_adam_dense_bias_mIdentity_60:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_60n
Identity_61IdentityRestoreV2:tensors:61"/device:CPU:0*
T0*
_output_shapes
:2
Identity_61?
AssignVariableOp_61AssignVariableOp)assignvariableop_61_adam_dense_1_kernel_mIdentity_61:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_61n
Identity_62IdentityRestoreV2:tensors:62"/device:CPU:0*
T0*
_output_shapes
:2
Identity_62?
AssignVariableOp_62AssignVariableOp'assignvariableop_62_adam_dense_1_bias_mIdentity_62:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_62n
Identity_63IdentityRestoreV2:tensors:63"/device:CPU:0*
T0*
_output_shapes
:2
Identity_63?
AssignVariableOp_63AssignVariableOp)assignvariableop_63_adam_dense_2_kernel_vIdentity_63:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_63n
Identity_64IdentityRestoreV2:tensors:64"/device:CPU:0*
T0*
_output_shapes
:2
Identity_64?
AssignVariableOp_64AssignVariableOp'assignvariableop_64_adam_dense_2_bias_vIdentity_64:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_64n
Identity_65IdentityRestoreV2:tensors:65"/device:CPU:0*
T0*
_output_shapes
:2
Identity_65?
AssignVariableOp_65AssignVariableOp(assignvariableop_65_adam_conv3d_kernel_vIdentity_65:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_65n
Identity_66IdentityRestoreV2:tensors:66"/device:CPU:0*
T0*
_output_shapes
:2
Identity_66?
AssignVariableOp_66AssignVariableOp&assignvariableop_66_adam_conv3d_bias_vIdentity_66:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_66n
Identity_67IdentityRestoreV2:tensors:67"/device:CPU:0*
T0*
_output_shapes
:2
Identity_67?
AssignVariableOp_67AssignVariableOp4assignvariableop_67_adam_batch_normalization_gamma_vIdentity_67:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_67n
Identity_68IdentityRestoreV2:tensors:68"/device:CPU:0*
T0*
_output_shapes
:2
Identity_68?
AssignVariableOp_68AssignVariableOp3assignvariableop_68_adam_batch_normalization_beta_vIdentity_68:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_68n
Identity_69IdentityRestoreV2:tensors:69"/device:CPU:0*
T0*
_output_shapes
:2
Identity_69?
AssignVariableOp_69AssignVariableOp*assignvariableop_69_adam_conv3d_1_kernel_vIdentity_69:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_69n
Identity_70IdentityRestoreV2:tensors:70"/device:CPU:0*
T0*
_output_shapes
:2
Identity_70?
AssignVariableOp_70AssignVariableOp(assignvariableop_70_adam_conv3d_1_bias_vIdentity_70:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_70n
Identity_71IdentityRestoreV2:tensors:71"/device:CPU:0*
T0*
_output_shapes
:2
Identity_71?
AssignVariableOp_71AssignVariableOp6assignvariableop_71_adam_batch_normalization_1_gamma_vIdentity_71:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_71n
Identity_72IdentityRestoreV2:tensors:72"/device:CPU:0*
T0*
_output_shapes
:2
Identity_72?
AssignVariableOp_72AssignVariableOp5assignvariableop_72_adam_batch_normalization_1_beta_vIdentity_72:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_72n
Identity_73IdentityRestoreV2:tensors:73"/device:CPU:0*
T0*
_output_shapes
:2
Identity_73?
AssignVariableOp_73AssignVariableOp*assignvariableop_73_adam_conv3d_2_kernel_vIdentity_73:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_73n
Identity_74IdentityRestoreV2:tensors:74"/device:CPU:0*
T0*
_output_shapes
:2
Identity_74?
AssignVariableOp_74AssignVariableOp(assignvariableop_74_adam_conv3d_2_bias_vIdentity_74:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_74n
Identity_75IdentityRestoreV2:tensors:75"/device:CPU:0*
T0*
_output_shapes
:2
Identity_75?
AssignVariableOp_75AssignVariableOp6assignvariableop_75_adam_batch_normalization_2_gamma_vIdentity_75:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_75n
Identity_76IdentityRestoreV2:tensors:76"/device:CPU:0*
T0*
_output_shapes
:2
Identity_76?
AssignVariableOp_76AssignVariableOp5assignvariableop_76_adam_batch_normalization_2_beta_vIdentity_76:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_76n
Identity_77IdentityRestoreV2:tensors:77"/device:CPU:0*
T0*
_output_shapes
:2
Identity_77?
AssignVariableOp_77AssignVariableOp*assignvariableop_77_adam_conv3d_3_kernel_vIdentity_77:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_77n
Identity_78IdentityRestoreV2:tensors:78"/device:CPU:0*
T0*
_output_shapes
:2
Identity_78?
AssignVariableOp_78AssignVariableOp(assignvariableop_78_adam_conv3d_3_bias_vIdentity_78:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_78n
Identity_79IdentityRestoreV2:tensors:79"/device:CPU:0*
T0*
_output_shapes
:2
Identity_79?
AssignVariableOp_79AssignVariableOp6assignvariableop_79_adam_batch_normalization_3_gamma_vIdentity_79:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_79n
Identity_80IdentityRestoreV2:tensors:80"/device:CPU:0*
T0*
_output_shapes
:2
Identity_80?
AssignVariableOp_80AssignVariableOp5assignvariableop_80_adam_batch_normalization_3_beta_vIdentity_80:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_80n
Identity_81IdentityRestoreV2:tensors:81"/device:CPU:0*
T0*
_output_shapes
:2
Identity_81?
AssignVariableOp_81AssignVariableOp(assignvariableop_81_adam_layer1_kernel_vIdentity_81:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_81n
Identity_82IdentityRestoreV2:tensors:82"/device:CPU:0*
T0*
_output_shapes
:2
Identity_82?
AssignVariableOp_82AssignVariableOp&assignvariableop_82_adam_layer1_bias_vIdentity_82:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_82n
Identity_83IdentityRestoreV2:tensors:83"/device:CPU:0*
T0*
_output_shapes
:2
Identity_83?
AssignVariableOp_83AssignVariableOp'assignvariableop_83_adam_dense_kernel_vIdentity_83:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_83n
Identity_84IdentityRestoreV2:tensors:84"/device:CPU:0*
T0*
_output_shapes
:2
Identity_84?
AssignVariableOp_84AssignVariableOp%assignvariableop_84_adam_dense_bias_vIdentity_84:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_84n
Identity_85IdentityRestoreV2:tensors:85"/device:CPU:0*
T0*
_output_shapes
:2
Identity_85?
AssignVariableOp_85AssignVariableOp)assignvariableop_85_adam_dense_1_kernel_vIdentity_85:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_85n
Identity_86IdentityRestoreV2:tensors:86"/device:CPU:0*
T0*
_output_shapes
:2
Identity_86?
AssignVariableOp_86AssignVariableOp'assignvariableop_86_adam_dense_1_bias_vIdentity_86:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_869
NoOpNoOp"/device:CPU:0*
_output_shapes
 2
NoOp?
Identity_87Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_49^AssignVariableOp_5^AssignVariableOp_50^AssignVariableOp_51^AssignVariableOp_52^AssignVariableOp_53^AssignVariableOp_54^AssignVariableOp_55^AssignVariableOp_56^AssignVariableOp_57^AssignVariableOp_58^AssignVariableOp_59^AssignVariableOp_6^AssignVariableOp_60^AssignVariableOp_61^AssignVariableOp_62^AssignVariableOp_63^AssignVariableOp_64^AssignVariableOp_65^AssignVariableOp_66^AssignVariableOp_67^AssignVariableOp_68^AssignVariableOp_69^AssignVariableOp_7^AssignVariableOp_70^AssignVariableOp_71^AssignVariableOp_72^AssignVariableOp_73^AssignVariableOp_74^AssignVariableOp_75^AssignVariableOp_76^AssignVariableOp_77^AssignVariableOp_78^AssignVariableOp_79^AssignVariableOp_8^AssignVariableOp_80^AssignVariableOp_81^AssignVariableOp_82^AssignVariableOp_83^AssignVariableOp_84^AssignVariableOp_85^AssignVariableOp_86^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: 2
Identity_87?
Identity_88IdentityIdentity_87:output:0^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_49^AssignVariableOp_5^AssignVariableOp_50^AssignVariableOp_51^AssignVariableOp_52^AssignVariableOp_53^AssignVariableOp_54^AssignVariableOp_55^AssignVariableOp_56^AssignVariableOp_57^AssignVariableOp_58^AssignVariableOp_59^AssignVariableOp_6^AssignVariableOp_60^AssignVariableOp_61^AssignVariableOp_62^AssignVariableOp_63^AssignVariableOp_64^AssignVariableOp_65^AssignVariableOp_66^AssignVariableOp_67^AssignVariableOp_68^AssignVariableOp_69^AssignVariableOp_7^AssignVariableOp_70^AssignVariableOp_71^AssignVariableOp_72^AssignVariableOp_73^AssignVariableOp_74^AssignVariableOp_75^AssignVariableOp_76^AssignVariableOp_77^AssignVariableOp_78^AssignVariableOp_79^AssignVariableOp_8^AssignVariableOp_80^AssignVariableOp_81^AssignVariableOp_82^AssignVariableOp_83^AssignVariableOp_84^AssignVariableOp_85^AssignVariableOp_86^AssignVariableOp_9*
T0*
_output_shapes
: 2
Identity_88"#
identity_88Identity_88:output:0*?
_input_shapes?
?: :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::2$
AssignVariableOpAssignVariableOp2(
AssignVariableOp_1AssignVariableOp_12*
AssignVariableOp_10AssignVariableOp_102*
AssignVariableOp_11AssignVariableOp_112*
AssignVariableOp_12AssignVariableOp_122*
AssignVariableOp_13AssignVariableOp_132*
AssignVariableOp_14AssignVariableOp_142*
AssignVariableOp_15AssignVariableOp_152*
AssignVariableOp_16AssignVariableOp_162*
AssignVariableOp_17AssignVariableOp_172*
AssignVariableOp_18AssignVariableOp_182*
AssignVariableOp_19AssignVariableOp_192(
AssignVariableOp_2AssignVariableOp_22*
AssignVariableOp_20AssignVariableOp_202*
AssignVariableOp_21AssignVariableOp_212*
AssignVariableOp_22AssignVariableOp_222*
AssignVariableOp_23AssignVariableOp_232*
AssignVariableOp_24AssignVariableOp_242*
AssignVariableOp_25AssignVariableOp_252*
AssignVariableOp_26AssignVariableOp_262*
AssignVariableOp_27AssignVariableOp_272*
AssignVariableOp_28AssignVariableOp_282*
AssignVariableOp_29AssignVariableOp_292(
AssignVariableOp_3AssignVariableOp_32*
AssignVariableOp_30AssignVariableOp_302*
AssignVariableOp_31AssignVariableOp_312*
AssignVariableOp_32AssignVariableOp_322*
AssignVariableOp_33AssignVariableOp_332*
AssignVariableOp_34AssignVariableOp_342*
AssignVariableOp_35AssignVariableOp_352*
AssignVariableOp_36AssignVariableOp_362*
AssignVariableOp_37AssignVariableOp_372*
AssignVariableOp_38AssignVariableOp_382*
AssignVariableOp_39AssignVariableOp_392(
AssignVariableOp_4AssignVariableOp_42*
AssignVariableOp_40AssignVariableOp_402*
AssignVariableOp_41AssignVariableOp_412*
AssignVariableOp_42AssignVariableOp_422*
AssignVariableOp_43AssignVariableOp_432*
AssignVariableOp_44AssignVariableOp_442*
AssignVariableOp_45AssignVariableOp_452*
AssignVariableOp_46AssignVariableOp_462*
AssignVariableOp_47AssignVariableOp_472*
AssignVariableOp_48AssignVariableOp_482*
AssignVariableOp_49AssignVariableOp_492(
AssignVariableOp_5AssignVariableOp_52*
AssignVariableOp_50AssignVariableOp_502*
AssignVariableOp_51AssignVariableOp_512*
AssignVariableOp_52AssignVariableOp_522*
AssignVariableOp_53AssignVariableOp_532*
AssignVariableOp_54AssignVariableOp_542*
AssignVariableOp_55AssignVariableOp_552*
AssignVariableOp_56AssignVariableOp_562*
AssignVariableOp_57AssignVariableOp_572*
AssignVariableOp_58AssignVariableOp_582*
AssignVariableOp_59AssignVariableOp_592(
AssignVariableOp_6AssignVariableOp_62*
AssignVariableOp_60AssignVariableOp_602*
AssignVariableOp_61AssignVariableOp_612*
AssignVariableOp_62AssignVariableOp_622*
AssignVariableOp_63AssignVariableOp_632*
AssignVariableOp_64AssignVariableOp_642*
AssignVariableOp_65AssignVariableOp_652*
AssignVariableOp_66AssignVariableOp_662*
AssignVariableOp_67AssignVariableOp_672*
AssignVariableOp_68AssignVariableOp_682*
AssignVariableOp_69AssignVariableOp_692(
AssignVariableOp_7AssignVariableOp_72*
AssignVariableOp_70AssignVariableOp_702*
AssignVariableOp_71AssignVariableOp_712*
AssignVariableOp_72AssignVariableOp_722*
AssignVariableOp_73AssignVariableOp_732*
AssignVariableOp_74AssignVariableOp_742*
AssignVariableOp_75AssignVariableOp_752*
AssignVariableOp_76AssignVariableOp_762*
AssignVariableOp_77AssignVariableOp_772*
AssignVariableOp_78AssignVariableOp_782*
AssignVariableOp_79AssignVariableOp_792(
AssignVariableOp_8AssignVariableOp_82*
AssignVariableOp_80AssignVariableOp_802*
AssignVariableOp_81AssignVariableOp_812*
AssignVariableOp_82AssignVariableOp_822*
AssignVariableOp_83AssignVariableOp_832*
AssignVariableOp_84AssignVariableOp_842*
AssignVariableOp_85AssignVariableOp_852*
AssignVariableOp_86AssignVariableOp_862(
AssignVariableOp_9AssignVariableOp_9:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
?
?
7__inference_batch_normalization_1_layer_call_fn_1227195

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *N
_output_shapes<
::8????????????????????????????????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_12224612
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????::::22
StatefulPartitionedCallStatefulPartitionedCall:v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
?
?
,__inference_sequential_layer_call_fn_1226304

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16

unknown_17

unknown_18

unknown_19

unknown_20

unknown_21

unknown_22

unknown_23

unknown_24

unknown_25

unknown_26

unknown_27

unknown_28
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19
unknown_20
unknown_21
unknown_22
unknown_23
unknown_24
unknown_25
unknown_26
unknown_27
unknown_28**
Tin#
!2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
*8
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *P
fKRI
G__inference_sequential_layer_call_and_return_conditional_losses_12240422
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????
2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????::::::::::::::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:\ X
4
_output_shapes"
 :??????????
 
_user_specified_nameinputs
?*
?
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_1227415

inputs
assignmovingavg_1227390
assignmovingavg_1_1227396)
%batchnorm_mul_readvariableop_resource%
!batchnorm_readvariableop_resource
identity??#AssignMovingAvg/AssignSubVariableOp?%AssignMovingAvg_1/AssignSubVariableOp?
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2 
moments/mean/reduction_indices?
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/mean?
moments/StopGradientStopGradientmoments/mean:output:0*
T0**
_output_shapes
:2
moments/StopGradient?
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*3
_output_shapes!
:?????????2
moments/SquaredDifference?
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2$
"moments/variance/reduction_indices?
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/variance?
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze?
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze_1?
AssignMovingAvg/decayConst**
_class 
loc:@AssignMovingAvg/1227390*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg/decay?
AssignMovingAvg/ReadVariableOpReadVariableOpassignmovingavg_1227390*
_output_shapes
:*
dtype02 
AssignMovingAvg/ReadVariableOp?
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0**
_class 
loc:@AssignMovingAvg/1227390*
_output_shapes
:2
AssignMovingAvg/sub?
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0**
_class 
loc:@AssignMovingAvg/1227390*
_output_shapes
:2
AssignMovingAvg/mul?
#AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1227390AssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp**
_class 
loc:@AssignMovingAvg/1227390*
_output_shapes
 *
dtype02%
#AssignMovingAvg/AssignSubVariableOp?
AssignMovingAvg_1/decayConst*,
_class"
 loc:@AssignMovingAvg_1/1227396*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg_1/decay?
 AssignMovingAvg_1/ReadVariableOpReadVariableOpassignmovingavg_1_1227396*
_output_shapes
:*
dtype02"
 AssignMovingAvg_1/ReadVariableOp?
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1227396*
_output_shapes
:2
AssignMovingAvg_1/sub?
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1227396*
_output_shapes
:2
AssignMovingAvg_1/mul?
%AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1_1227396AssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*,
_class"
 loc:@AssignMovingAvg_1/1227396*
_output_shapes
 *
dtype02'
%AssignMovingAvg_1/AssignSubVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/mul_1{
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp?
batchnorm/subSub batchnorm/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0$^AssignMovingAvg/AssignSubVariableOp&^AssignMovingAvg_1/AssignSubVariableOp*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????::::2J
#AssignMovingAvg/AssignSubVariableOp#AssignMovingAvg/AssignSubVariableOp2N
%AssignMovingAvg_1/AssignSubVariableOp%AssignMovingAvg_1/AssignSubVariableOp:[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?+
?
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_1227333

inputs
assignmovingavg_1227308
assignmovingavg_1_1227314)
%batchnorm_mul_readvariableop_resource%
!batchnorm_readvariableop_resource
identity??#AssignMovingAvg/AssignSubVariableOp?%AssignMovingAvg_1/AssignSubVariableOp?
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2 
moments/mean/reduction_indices?
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/mean?
moments/StopGradientStopGradientmoments/mean:output:0*
T0**
_output_shapes
:2
moments/StopGradient?
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
moments/SquaredDifference?
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2$
"moments/variance/reduction_indices?
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/variance?
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze?
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze_1?
AssignMovingAvg/decayConst**
_class 
loc:@AssignMovingAvg/1227308*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg/decay?
AssignMovingAvg/ReadVariableOpReadVariableOpassignmovingavg_1227308*
_output_shapes
:*
dtype02 
AssignMovingAvg/ReadVariableOp?
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0**
_class 
loc:@AssignMovingAvg/1227308*
_output_shapes
:2
AssignMovingAvg/sub?
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0**
_class 
loc:@AssignMovingAvg/1227308*
_output_shapes
:2
AssignMovingAvg/mul?
#AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1227308AssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp**
_class 
loc:@AssignMovingAvg/1227308*
_output_shapes
 *
dtype02%
#AssignMovingAvg/AssignSubVariableOp?
AssignMovingAvg_1/decayConst*,
_class"
 loc:@AssignMovingAvg_1/1227314*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg_1/decay?
 AssignMovingAvg_1/ReadVariableOpReadVariableOpassignmovingavg_1_1227314*
_output_shapes
:*
dtype02"
 AssignMovingAvg_1/ReadVariableOp?
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1227314*
_output_shapes
:2
AssignMovingAvg_1/sub?
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1227314*
_output_shapes
:2
AssignMovingAvg_1/mul?
%AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1_1227314AssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*,
_class"
 loc:@AssignMovingAvg_1/1227314*
_output_shapes
 *
dtype02'
%AssignMovingAvg_1/AssignSubVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/mul_1{
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp?
batchnorm/subSub batchnorm/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0$^AssignMovingAvg/AssignSubVariableOp&^AssignMovingAvg_1/AssignSubVariableOp*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????::::2J
#AssignMovingAvg/AssignSubVariableOp#AssignMovingAvg/AssignSubVariableOp2N
%AssignMovingAvg_1/AssignSubVariableOp%AssignMovingAvg_1/AssignSubVariableOp:v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_1222967

inputs%
!batchnorm_readvariableop_resource)
%batchnorm_mul_readvariableop_resource'
#batchnorm_readvariableop_1_resource'
#batchnorm_readvariableop_2_resource
identity??
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/mul_1?
batchnorm/ReadVariableOp_1ReadVariableOp#batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_1?
batchnorm/mul_2Mul"batchnorm/ReadVariableOp_1:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOp_2ReadVariableOp#batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_2?
batchnorm/subSub"batchnorm/ReadVariableOp_2:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/add_1s
IdentityIdentitybatchnorm/add_1:z:0*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????:::::[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?C
?	
I__inference_functional_1_layer_call_and_return_conditional_losses_1223451

inputs
conv3d_1223386
conv3d_1223388
batch_normalization_1223391
batch_normalization_1223393
batch_normalization_1223395
batch_normalization_1223397
conv3d_1_1223400
conv3d_1_1223402!
batch_normalization_1_1223405!
batch_normalization_1_1223407!
batch_normalization_1_1223409!
batch_normalization_1_1223411
conv3d_2_1223414
conv3d_2_1223416!
batch_normalization_2_1223419!
batch_normalization_2_1223421!
batch_normalization_2_1223423!
batch_normalization_2_1223425
conv3d_3_1223428
conv3d_3_1223430!
batch_normalization_3_1223433!
batch_normalization_3_1223435!
batch_normalization_3_1223437!
batch_normalization_3_1223439
layer1_1223445
layer1_1223447
identity??+batch_normalization/StatefulPartitionedCall?-batch_normalization_1/StatefulPartitionedCall?-batch_normalization_2/StatefulPartitionedCall?-batch_normalization_3/StatefulPartitionedCall?conv3d/StatefulPartitionedCall? conv3d_1/StatefulPartitionedCall? conv3d_2/StatefulPartitionedCall? conv3d_3/StatefulPartitionedCall?dropout/StatefulPartitionedCall?layer1/StatefulPartitionedCall?
conv3d/StatefulPartitionedCallStatefulPartitionedCallinputsconv3d_1223386conv3d_1223388*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *L
fGRE
C__inference_conv3d_layer_call_and_return_conditional_losses_12227782 
conv3d/StatefulPartitionedCall?
+batch_normalization/StatefulPartitionedCallStatefulPartitionedCall'conv3d/StatefulPartitionedCall:output:0batch_normalization_1223391batch_normalization_1223393batch_normalization_1223395batch_normalization_1223397*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *Y
fTRR
P__inference_batch_normalization_layer_call_and_return_conditional_losses_12228292-
+batch_normalization/StatefulPartitionedCall?
 conv3d_1/StatefulPartitionedCallStatefulPartitionedCall4batch_normalization/StatefulPartitionedCall:output:0conv3d_1_1223400conv3d_1_1223402*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv3d_1_layer_call_and_return_conditional_losses_12228962"
 conv3d_1/StatefulPartitionedCall?
-batch_normalization_1/StatefulPartitionedCallStatefulPartitionedCall)conv3d_1/StatefulPartitionedCall:output:0batch_normalization_1_1223405batch_normalization_1_1223407batch_normalization_1_1223409batch_normalization_1_1223411*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_12229472/
-batch_normalization_1/StatefulPartitionedCall?
 conv3d_2/StatefulPartitionedCallStatefulPartitionedCall6batch_normalization_1/StatefulPartitionedCall:output:0conv3d_2_1223414conv3d_2_1223416*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv3d_2_layer_call_and_return_conditional_losses_12230142"
 conv3d_2/StatefulPartitionedCall?
-batch_normalization_2/StatefulPartitionedCallStatefulPartitionedCall)conv3d_2/StatefulPartitionedCall:output:0batch_normalization_2_1223419batch_normalization_2_1223421batch_normalization_2_1223423batch_normalization_2_1223425*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_12230652/
-batch_normalization_2/StatefulPartitionedCall?
 conv3d_3/StatefulPartitionedCallStatefulPartitionedCall6batch_normalization_2/StatefulPartitionedCall:output:0conv3d_3_1223428conv3d_3_1223430*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv3d_3_layer_call_and_return_conditional_losses_12231322"
 conv3d_3/StatefulPartitionedCall?
-batch_normalization_3/StatefulPartitionedCallStatefulPartitionedCall)conv3d_3/StatefulPartitionedCall:output:0batch_normalization_3_1223433batch_normalization_3_1223435batch_normalization_3_1223437batch_normalization_3_1223439*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_12231832/
-batch_normalization_3/StatefulPartitionedCall?
!average_pooling3d/PartitionedCallPartitionedCall6batch_normalization_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *W
fRRP
N__inference_average_pooling3d_layer_call_and_return_conditional_losses_12227582#
!average_pooling3d/PartitionedCall?
flatten/PartitionedCallPartitionedCall*average_pooling3d/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????
* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_flatten_layer_call_and_return_conditional_losses_12232462
flatten/PartitionedCall?
dropout/StatefulPartitionedCallStatefulPartitionedCall flatten/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????
* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dropout_layer_call_and_return_conditional_losses_12232662!
dropout/StatefulPartitionedCall?
layer1/StatefulPartitionedCallStatefulPartitionedCall(dropout/StatefulPartitionedCall:output:0layer1_1223445layer1_1223447*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *L
fGRE
C__inference_layer1_layer_call_and_return_conditional_losses_12232952 
layer1/StatefulPartitionedCall?
IdentityIdentity'layer1/StatefulPartitionedCall:output:0,^batch_normalization/StatefulPartitionedCall.^batch_normalization_1/StatefulPartitionedCall.^batch_normalization_2/StatefulPartitionedCall.^batch_normalization_3/StatefulPartitionedCall^conv3d/StatefulPartitionedCall!^conv3d_1/StatefulPartitionedCall!^conv3d_2/StatefulPartitionedCall!^conv3d_3/StatefulPartitionedCall ^dropout/StatefulPartitionedCall^layer1/StatefulPartitionedCall*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????::::::::::::::::::::::::::2Z
+batch_normalization/StatefulPartitionedCall+batch_normalization/StatefulPartitionedCall2^
-batch_normalization_1/StatefulPartitionedCall-batch_normalization_1/StatefulPartitionedCall2^
-batch_normalization_2/StatefulPartitionedCall-batch_normalization_2/StatefulPartitionedCall2^
-batch_normalization_3/StatefulPartitionedCall-batch_normalization_3/StatefulPartitionedCall2@
conv3d/StatefulPartitionedCallconv3d/StatefulPartitionedCall2D
 conv3d_1/StatefulPartitionedCall conv3d_1/StatefulPartitionedCall2D
 conv3d_2/StatefulPartitionedCall conv3d_2/StatefulPartitionedCall2D
 conv3d_3/StatefulPartitionedCall conv3d_3/StatefulPartitionedCall2B
dropout/StatefulPartitionedCalldropout/StatefulPartitionedCall2@
layer1/StatefulPartitionedCalllayer1/StatefulPartitionedCall:\ X
4
_output_shapes"
 :??????????
 
_user_specified_nameinputs
?
?
5__inference_batch_normalization_layer_call_fn_1227093

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *Y
fTRR
P__inference_batch_normalization_layer_call_and_return_conditional_losses_12228492
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????::::22
StatefulPartitionedCallStatefulPartitionedCall:[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?
`
D__inference_flatten_layer_call_and_return_conditional_losses_1223246

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"????   2
Consth
ReshapeReshapeinputsConst:output:0*
T0*(
_output_shapes
:??????????
2	
Reshapee
IdentityIdentityReshape:output:0*
T0*(
_output_shapes
:??????????
2

Identity"
identityIdentity:output:0*2
_input_shapes!
:?????????:[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
??
?
I__inference_functional_1_layer_call_and_return_conditional_losses_1226604

inputs)
%conv3d_conv3d_readvariableop_resource*
&conv3d_biasadd_readvariableop_resource/
+batch_normalization_assignmovingavg_12264441
-batch_normalization_assignmovingavg_1_1226450=
9batch_normalization_batchnorm_mul_readvariableop_resource9
5batch_normalization_batchnorm_readvariableop_resource+
'conv3d_1_conv3d_readvariableop_resource,
(conv3d_1_biasadd_readvariableop_resource1
-batch_normalization_1_assignmovingavg_12264833
/batch_normalization_1_assignmovingavg_1_1226489?
;batch_normalization_1_batchnorm_mul_readvariableop_resource;
7batch_normalization_1_batchnorm_readvariableop_resource+
'conv3d_2_conv3d_readvariableop_resource,
(conv3d_2_biasadd_readvariableop_resource1
-batch_normalization_2_assignmovingavg_12265223
/batch_normalization_2_assignmovingavg_1_1226528?
;batch_normalization_2_batchnorm_mul_readvariableop_resource;
7batch_normalization_2_batchnorm_readvariableop_resource+
'conv3d_3_conv3d_readvariableop_resource,
(conv3d_3_biasadd_readvariableop_resource1
-batch_normalization_3_assignmovingavg_12265613
/batch_normalization_3_assignmovingavg_1_1226567?
;batch_normalization_3_batchnorm_mul_readvariableop_resource;
7batch_normalization_3_batchnorm_readvariableop_resource)
%layer1_matmul_readvariableop_resource*
&layer1_biasadd_readvariableop_resource
identity??7batch_normalization/AssignMovingAvg/AssignSubVariableOp?9batch_normalization/AssignMovingAvg_1/AssignSubVariableOp?9batch_normalization_1/AssignMovingAvg/AssignSubVariableOp?;batch_normalization_1/AssignMovingAvg_1/AssignSubVariableOp?9batch_normalization_2/AssignMovingAvg/AssignSubVariableOp?;batch_normalization_2/AssignMovingAvg_1/AssignSubVariableOp?9batch_normalization_3/AssignMovingAvg/AssignSubVariableOp?;batch_normalization_3/AssignMovingAvg_1/AssignSubVariableOp?
conv3d/Conv3D/ReadVariableOpReadVariableOp%conv3d_conv3d_readvariableop_resource*+
_output_shapes
:?*
dtype02
conv3d/Conv3D/ReadVariableOp?
conv3d/Conv3DConv3Dinputs$conv3d/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
conv3d/Conv3D?
conv3d/BiasAdd/ReadVariableOpReadVariableOp&conv3d_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02
conv3d/BiasAdd/ReadVariableOp?
conv3d/BiasAddBiasAddconv3d/Conv3D:output:0%conv3d/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2
conv3d/BiasAdd?
2batch_normalization/moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             24
2batch_normalization/moments/mean/reduction_indices?
 batch_normalization/moments/meanMeanconv3d/BiasAdd:output:0;batch_normalization/moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2"
 batch_normalization/moments/mean?
(batch_normalization/moments/StopGradientStopGradient)batch_normalization/moments/mean:output:0*
T0**
_output_shapes
:2*
(batch_normalization/moments/StopGradient?
-batch_normalization/moments/SquaredDifferenceSquaredDifferenceconv3d/BiasAdd:output:01batch_normalization/moments/StopGradient:output:0*
T0*3
_output_shapes!
:?????????2/
-batch_normalization/moments/SquaredDifference?
6batch_normalization/moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             28
6batch_normalization/moments/variance/reduction_indices?
$batch_normalization/moments/varianceMean1batch_normalization/moments/SquaredDifference:z:0?batch_normalization/moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2&
$batch_normalization/moments/variance?
#batch_normalization/moments/SqueezeSqueeze)batch_normalization/moments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2%
#batch_normalization/moments/Squeeze?
%batch_normalization/moments/Squeeze_1Squeeze-batch_normalization/moments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2'
%batch_normalization/moments/Squeeze_1?
)batch_normalization/AssignMovingAvg/decayConst*>
_class4
20loc:@batch_normalization/AssignMovingAvg/1226444*
_output_shapes
: *
dtype0*
valueB
 *
?#<2+
)batch_normalization/AssignMovingAvg/decay?
2batch_normalization/AssignMovingAvg/ReadVariableOpReadVariableOp+batch_normalization_assignmovingavg_1226444*
_output_shapes
:*
dtype024
2batch_normalization/AssignMovingAvg/ReadVariableOp?
'batch_normalization/AssignMovingAvg/subSub:batch_normalization/AssignMovingAvg/ReadVariableOp:value:0,batch_normalization/moments/Squeeze:output:0*
T0*>
_class4
20loc:@batch_normalization/AssignMovingAvg/1226444*
_output_shapes
:2)
'batch_normalization/AssignMovingAvg/sub?
'batch_normalization/AssignMovingAvg/mulMul+batch_normalization/AssignMovingAvg/sub:z:02batch_normalization/AssignMovingAvg/decay:output:0*
T0*>
_class4
20loc:@batch_normalization/AssignMovingAvg/1226444*
_output_shapes
:2)
'batch_normalization/AssignMovingAvg/mul?
7batch_normalization/AssignMovingAvg/AssignSubVariableOpAssignSubVariableOp+batch_normalization_assignmovingavg_1226444+batch_normalization/AssignMovingAvg/mul:z:03^batch_normalization/AssignMovingAvg/ReadVariableOp*>
_class4
20loc:@batch_normalization/AssignMovingAvg/1226444*
_output_shapes
 *
dtype029
7batch_normalization/AssignMovingAvg/AssignSubVariableOp?
+batch_normalization/AssignMovingAvg_1/decayConst*@
_class6
42loc:@batch_normalization/AssignMovingAvg_1/1226450*
_output_shapes
: *
dtype0*
valueB
 *
?#<2-
+batch_normalization/AssignMovingAvg_1/decay?
4batch_normalization/AssignMovingAvg_1/ReadVariableOpReadVariableOp-batch_normalization_assignmovingavg_1_1226450*
_output_shapes
:*
dtype026
4batch_normalization/AssignMovingAvg_1/ReadVariableOp?
)batch_normalization/AssignMovingAvg_1/subSub<batch_normalization/AssignMovingAvg_1/ReadVariableOp:value:0.batch_normalization/moments/Squeeze_1:output:0*
T0*@
_class6
42loc:@batch_normalization/AssignMovingAvg_1/1226450*
_output_shapes
:2+
)batch_normalization/AssignMovingAvg_1/sub?
)batch_normalization/AssignMovingAvg_1/mulMul-batch_normalization/AssignMovingAvg_1/sub:z:04batch_normalization/AssignMovingAvg_1/decay:output:0*
T0*@
_class6
42loc:@batch_normalization/AssignMovingAvg_1/1226450*
_output_shapes
:2+
)batch_normalization/AssignMovingAvg_1/mul?
9batch_normalization/AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOp-batch_normalization_assignmovingavg_1_1226450-batch_normalization/AssignMovingAvg_1/mul:z:05^batch_normalization/AssignMovingAvg_1/ReadVariableOp*@
_class6
42loc:@batch_normalization/AssignMovingAvg_1/1226450*
_output_shapes
 *
dtype02;
9batch_normalization/AssignMovingAvg_1/AssignSubVariableOp?
#batch_normalization/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2%
#batch_normalization/batchnorm/add/y?
!batch_normalization/batchnorm/addAddV2.batch_normalization/moments/Squeeze_1:output:0,batch_normalization/batchnorm/add/y:output:0*
T0*
_output_shapes
:2#
!batch_normalization/batchnorm/add?
#batch_normalization/batchnorm/RsqrtRsqrt%batch_normalization/batchnorm/add:z:0*
T0*
_output_shapes
:2%
#batch_normalization/batchnorm/Rsqrt?
0batch_normalization/batchnorm/mul/ReadVariableOpReadVariableOp9batch_normalization_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype022
0batch_normalization/batchnorm/mul/ReadVariableOp?
!batch_normalization/batchnorm/mulMul'batch_normalization/batchnorm/Rsqrt:y:08batch_normalization/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2#
!batch_normalization/batchnorm/mul?
#batch_normalization/batchnorm/mul_1Mulconv3d/BiasAdd:output:0%batch_normalization/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2%
#batch_normalization/batchnorm/mul_1?
#batch_normalization/batchnorm/mul_2Mul,batch_normalization/moments/Squeeze:output:0%batch_normalization/batchnorm/mul:z:0*
T0*
_output_shapes
:2%
#batch_normalization/batchnorm/mul_2?
,batch_normalization/batchnorm/ReadVariableOpReadVariableOp5batch_normalization_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02.
,batch_normalization/batchnorm/ReadVariableOp?
!batch_normalization/batchnorm/subSub4batch_normalization/batchnorm/ReadVariableOp:value:0'batch_normalization/batchnorm/mul_2:z:0*
T0*
_output_shapes
:2#
!batch_normalization/batchnorm/sub?
#batch_normalization/batchnorm/add_1AddV2'batch_normalization/batchnorm/mul_1:z:0%batch_normalization/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2%
#batch_normalization/batchnorm/add_1?
conv3d_1/Conv3D/ReadVariableOpReadVariableOp'conv3d_1_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02 
conv3d_1/Conv3D/ReadVariableOp?
conv3d_1/Conv3DConv3D'batch_normalization/batchnorm/add_1:z:0&conv3d_1/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
conv3d_1/Conv3D?
conv3d_1/BiasAdd/ReadVariableOpReadVariableOp(conv3d_1_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02!
conv3d_1/BiasAdd/ReadVariableOp?
conv3d_1/BiasAddBiasAddconv3d_1/Conv3D:output:0'conv3d_1/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2
conv3d_1/BiasAdd|
conv3d_1/EluEluconv3d_1/BiasAdd:output:0*
T0*3
_output_shapes!
:?????????2
conv3d_1/Elu?
4batch_normalization_1/moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             26
4batch_normalization_1/moments/mean/reduction_indices?
"batch_normalization_1/moments/meanMeanconv3d_1/Elu:activations:0=batch_normalization_1/moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2$
"batch_normalization_1/moments/mean?
*batch_normalization_1/moments/StopGradientStopGradient+batch_normalization_1/moments/mean:output:0*
T0**
_output_shapes
:2,
*batch_normalization_1/moments/StopGradient?
/batch_normalization_1/moments/SquaredDifferenceSquaredDifferenceconv3d_1/Elu:activations:03batch_normalization_1/moments/StopGradient:output:0*
T0*3
_output_shapes!
:?????????21
/batch_normalization_1/moments/SquaredDifference?
8batch_normalization_1/moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2:
8batch_normalization_1/moments/variance/reduction_indices?
&batch_normalization_1/moments/varianceMean3batch_normalization_1/moments/SquaredDifference:z:0Abatch_normalization_1/moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2(
&batch_normalization_1/moments/variance?
%batch_normalization_1/moments/SqueezeSqueeze+batch_normalization_1/moments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2'
%batch_normalization_1/moments/Squeeze?
'batch_normalization_1/moments/Squeeze_1Squeeze/batch_normalization_1/moments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2)
'batch_normalization_1/moments/Squeeze_1?
+batch_normalization_1/AssignMovingAvg/decayConst*@
_class6
42loc:@batch_normalization_1/AssignMovingAvg/1226483*
_output_shapes
: *
dtype0*
valueB
 *
?#<2-
+batch_normalization_1/AssignMovingAvg/decay?
4batch_normalization_1/AssignMovingAvg/ReadVariableOpReadVariableOp-batch_normalization_1_assignmovingavg_1226483*
_output_shapes
:*
dtype026
4batch_normalization_1/AssignMovingAvg/ReadVariableOp?
)batch_normalization_1/AssignMovingAvg/subSub<batch_normalization_1/AssignMovingAvg/ReadVariableOp:value:0.batch_normalization_1/moments/Squeeze:output:0*
T0*@
_class6
42loc:@batch_normalization_1/AssignMovingAvg/1226483*
_output_shapes
:2+
)batch_normalization_1/AssignMovingAvg/sub?
)batch_normalization_1/AssignMovingAvg/mulMul-batch_normalization_1/AssignMovingAvg/sub:z:04batch_normalization_1/AssignMovingAvg/decay:output:0*
T0*@
_class6
42loc:@batch_normalization_1/AssignMovingAvg/1226483*
_output_shapes
:2+
)batch_normalization_1/AssignMovingAvg/mul?
9batch_normalization_1/AssignMovingAvg/AssignSubVariableOpAssignSubVariableOp-batch_normalization_1_assignmovingavg_1226483-batch_normalization_1/AssignMovingAvg/mul:z:05^batch_normalization_1/AssignMovingAvg/ReadVariableOp*@
_class6
42loc:@batch_normalization_1/AssignMovingAvg/1226483*
_output_shapes
 *
dtype02;
9batch_normalization_1/AssignMovingAvg/AssignSubVariableOp?
-batch_normalization_1/AssignMovingAvg_1/decayConst*B
_class8
64loc:@batch_normalization_1/AssignMovingAvg_1/1226489*
_output_shapes
: *
dtype0*
valueB
 *
?#<2/
-batch_normalization_1/AssignMovingAvg_1/decay?
6batch_normalization_1/AssignMovingAvg_1/ReadVariableOpReadVariableOp/batch_normalization_1_assignmovingavg_1_1226489*
_output_shapes
:*
dtype028
6batch_normalization_1/AssignMovingAvg_1/ReadVariableOp?
+batch_normalization_1/AssignMovingAvg_1/subSub>batch_normalization_1/AssignMovingAvg_1/ReadVariableOp:value:00batch_normalization_1/moments/Squeeze_1:output:0*
T0*B
_class8
64loc:@batch_normalization_1/AssignMovingAvg_1/1226489*
_output_shapes
:2-
+batch_normalization_1/AssignMovingAvg_1/sub?
+batch_normalization_1/AssignMovingAvg_1/mulMul/batch_normalization_1/AssignMovingAvg_1/sub:z:06batch_normalization_1/AssignMovingAvg_1/decay:output:0*
T0*B
_class8
64loc:@batch_normalization_1/AssignMovingAvg_1/1226489*
_output_shapes
:2-
+batch_normalization_1/AssignMovingAvg_1/mul?
;batch_normalization_1/AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOp/batch_normalization_1_assignmovingavg_1_1226489/batch_normalization_1/AssignMovingAvg_1/mul:z:07^batch_normalization_1/AssignMovingAvg_1/ReadVariableOp*B
_class8
64loc:@batch_normalization_1/AssignMovingAvg_1/1226489*
_output_shapes
 *
dtype02=
;batch_normalization_1/AssignMovingAvg_1/AssignSubVariableOp?
%batch_normalization_1/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2'
%batch_normalization_1/batchnorm/add/y?
#batch_normalization_1/batchnorm/addAddV20batch_normalization_1/moments/Squeeze_1:output:0.batch_normalization_1/batchnorm/add/y:output:0*
T0*
_output_shapes
:2%
#batch_normalization_1/batchnorm/add?
%batch_normalization_1/batchnorm/RsqrtRsqrt'batch_normalization_1/batchnorm/add:z:0*
T0*
_output_shapes
:2'
%batch_normalization_1/batchnorm/Rsqrt?
2batch_normalization_1/batchnorm/mul/ReadVariableOpReadVariableOp;batch_normalization_1_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype024
2batch_normalization_1/batchnorm/mul/ReadVariableOp?
#batch_normalization_1/batchnorm/mulMul)batch_normalization_1/batchnorm/Rsqrt:y:0:batch_normalization_1/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2%
#batch_normalization_1/batchnorm/mul?
%batch_normalization_1/batchnorm/mul_1Mulconv3d_1/Elu:activations:0'batch_normalization_1/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2'
%batch_normalization_1/batchnorm/mul_1?
%batch_normalization_1/batchnorm/mul_2Mul.batch_normalization_1/moments/Squeeze:output:0'batch_normalization_1/batchnorm/mul:z:0*
T0*
_output_shapes
:2'
%batch_normalization_1/batchnorm/mul_2?
.batch_normalization_1/batchnorm/ReadVariableOpReadVariableOp7batch_normalization_1_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype020
.batch_normalization_1/batchnorm/ReadVariableOp?
#batch_normalization_1/batchnorm/subSub6batch_normalization_1/batchnorm/ReadVariableOp:value:0)batch_normalization_1/batchnorm/mul_2:z:0*
T0*
_output_shapes
:2%
#batch_normalization_1/batchnorm/sub?
%batch_normalization_1/batchnorm/add_1AddV2)batch_normalization_1/batchnorm/mul_1:z:0'batch_normalization_1/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2'
%batch_normalization_1/batchnorm/add_1?
conv3d_2/Conv3D/ReadVariableOpReadVariableOp'conv3d_2_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02 
conv3d_2/Conv3D/ReadVariableOp?
conv3d_2/Conv3DConv3D)batch_normalization_1/batchnorm/add_1:z:0&conv3d_2/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
conv3d_2/Conv3D?
conv3d_2/BiasAdd/ReadVariableOpReadVariableOp(conv3d_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02!
conv3d_2/BiasAdd/ReadVariableOp?
conv3d_2/BiasAddBiasAddconv3d_2/Conv3D:output:0'conv3d_2/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2
conv3d_2/BiasAdd|
conv3d_2/EluEluconv3d_2/BiasAdd:output:0*
T0*3
_output_shapes!
:?????????2
conv3d_2/Elu?
4batch_normalization_2/moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             26
4batch_normalization_2/moments/mean/reduction_indices?
"batch_normalization_2/moments/meanMeanconv3d_2/Elu:activations:0=batch_normalization_2/moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2$
"batch_normalization_2/moments/mean?
*batch_normalization_2/moments/StopGradientStopGradient+batch_normalization_2/moments/mean:output:0*
T0**
_output_shapes
:2,
*batch_normalization_2/moments/StopGradient?
/batch_normalization_2/moments/SquaredDifferenceSquaredDifferenceconv3d_2/Elu:activations:03batch_normalization_2/moments/StopGradient:output:0*
T0*3
_output_shapes!
:?????????21
/batch_normalization_2/moments/SquaredDifference?
8batch_normalization_2/moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2:
8batch_normalization_2/moments/variance/reduction_indices?
&batch_normalization_2/moments/varianceMean3batch_normalization_2/moments/SquaredDifference:z:0Abatch_normalization_2/moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2(
&batch_normalization_2/moments/variance?
%batch_normalization_2/moments/SqueezeSqueeze+batch_normalization_2/moments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2'
%batch_normalization_2/moments/Squeeze?
'batch_normalization_2/moments/Squeeze_1Squeeze/batch_normalization_2/moments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2)
'batch_normalization_2/moments/Squeeze_1?
+batch_normalization_2/AssignMovingAvg/decayConst*@
_class6
42loc:@batch_normalization_2/AssignMovingAvg/1226522*
_output_shapes
: *
dtype0*
valueB
 *
?#<2-
+batch_normalization_2/AssignMovingAvg/decay?
4batch_normalization_2/AssignMovingAvg/ReadVariableOpReadVariableOp-batch_normalization_2_assignmovingavg_1226522*
_output_shapes
:*
dtype026
4batch_normalization_2/AssignMovingAvg/ReadVariableOp?
)batch_normalization_2/AssignMovingAvg/subSub<batch_normalization_2/AssignMovingAvg/ReadVariableOp:value:0.batch_normalization_2/moments/Squeeze:output:0*
T0*@
_class6
42loc:@batch_normalization_2/AssignMovingAvg/1226522*
_output_shapes
:2+
)batch_normalization_2/AssignMovingAvg/sub?
)batch_normalization_2/AssignMovingAvg/mulMul-batch_normalization_2/AssignMovingAvg/sub:z:04batch_normalization_2/AssignMovingAvg/decay:output:0*
T0*@
_class6
42loc:@batch_normalization_2/AssignMovingAvg/1226522*
_output_shapes
:2+
)batch_normalization_2/AssignMovingAvg/mul?
9batch_normalization_2/AssignMovingAvg/AssignSubVariableOpAssignSubVariableOp-batch_normalization_2_assignmovingavg_1226522-batch_normalization_2/AssignMovingAvg/mul:z:05^batch_normalization_2/AssignMovingAvg/ReadVariableOp*@
_class6
42loc:@batch_normalization_2/AssignMovingAvg/1226522*
_output_shapes
 *
dtype02;
9batch_normalization_2/AssignMovingAvg/AssignSubVariableOp?
-batch_normalization_2/AssignMovingAvg_1/decayConst*B
_class8
64loc:@batch_normalization_2/AssignMovingAvg_1/1226528*
_output_shapes
: *
dtype0*
valueB
 *
?#<2/
-batch_normalization_2/AssignMovingAvg_1/decay?
6batch_normalization_2/AssignMovingAvg_1/ReadVariableOpReadVariableOp/batch_normalization_2_assignmovingavg_1_1226528*
_output_shapes
:*
dtype028
6batch_normalization_2/AssignMovingAvg_1/ReadVariableOp?
+batch_normalization_2/AssignMovingAvg_1/subSub>batch_normalization_2/AssignMovingAvg_1/ReadVariableOp:value:00batch_normalization_2/moments/Squeeze_1:output:0*
T0*B
_class8
64loc:@batch_normalization_2/AssignMovingAvg_1/1226528*
_output_shapes
:2-
+batch_normalization_2/AssignMovingAvg_1/sub?
+batch_normalization_2/AssignMovingAvg_1/mulMul/batch_normalization_2/AssignMovingAvg_1/sub:z:06batch_normalization_2/AssignMovingAvg_1/decay:output:0*
T0*B
_class8
64loc:@batch_normalization_2/AssignMovingAvg_1/1226528*
_output_shapes
:2-
+batch_normalization_2/AssignMovingAvg_1/mul?
;batch_normalization_2/AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOp/batch_normalization_2_assignmovingavg_1_1226528/batch_normalization_2/AssignMovingAvg_1/mul:z:07^batch_normalization_2/AssignMovingAvg_1/ReadVariableOp*B
_class8
64loc:@batch_normalization_2/AssignMovingAvg_1/1226528*
_output_shapes
 *
dtype02=
;batch_normalization_2/AssignMovingAvg_1/AssignSubVariableOp?
%batch_normalization_2/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2'
%batch_normalization_2/batchnorm/add/y?
#batch_normalization_2/batchnorm/addAddV20batch_normalization_2/moments/Squeeze_1:output:0.batch_normalization_2/batchnorm/add/y:output:0*
T0*
_output_shapes
:2%
#batch_normalization_2/batchnorm/add?
%batch_normalization_2/batchnorm/RsqrtRsqrt'batch_normalization_2/batchnorm/add:z:0*
T0*
_output_shapes
:2'
%batch_normalization_2/batchnorm/Rsqrt?
2batch_normalization_2/batchnorm/mul/ReadVariableOpReadVariableOp;batch_normalization_2_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype024
2batch_normalization_2/batchnorm/mul/ReadVariableOp?
#batch_normalization_2/batchnorm/mulMul)batch_normalization_2/batchnorm/Rsqrt:y:0:batch_normalization_2/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2%
#batch_normalization_2/batchnorm/mul?
%batch_normalization_2/batchnorm/mul_1Mulconv3d_2/Elu:activations:0'batch_normalization_2/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2'
%batch_normalization_2/batchnorm/mul_1?
%batch_normalization_2/batchnorm/mul_2Mul.batch_normalization_2/moments/Squeeze:output:0'batch_normalization_2/batchnorm/mul:z:0*
T0*
_output_shapes
:2'
%batch_normalization_2/batchnorm/mul_2?
.batch_normalization_2/batchnorm/ReadVariableOpReadVariableOp7batch_normalization_2_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype020
.batch_normalization_2/batchnorm/ReadVariableOp?
#batch_normalization_2/batchnorm/subSub6batch_normalization_2/batchnorm/ReadVariableOp:value:0)batch_normalization_2/batchnorm/mul_2:z:0*
T0*
_output_shapes
:2%
#batch_normalization_2/batchnorm/sub?
%batch_normalization_2/batchnorm/add_1AddV2)batch_normalization_2/batchnorm/mul_1:z:0'batch_normalization_2/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2'
%batch_normalization_2/batchnorm/add_1?
conv3d_3/Conv3D/ReadVariableOpReadVariableOp'conv3d_3_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02 
conv3d_3/Conv3D/ReadVariableOp?
conv3d_3/Conv3DConv3D)batch_normalization_2/batchnorm/add_1:z:0&conv3d_3/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
conv3d_3/Conv3D?
conv3d_3/BiasAdd/ReadVariableOpReadVariableOp(conv3d_3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02!
conv3d_3/BiasAdd/ReadVariableOp?
conv3d_3/BiasAddBiasAddconv3d_3/Conv3D:output:0'conv3d_3/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2
conv3d_3/BiasAdd|
conv3d_3/EluEluconv3d_3/BiasAdd:output:0*
T0*3
_output_shapes!
:?????????2
conv3d_3/Elu?
4batch_normalization_3/moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             26
4batch_normalization_3/moments/mean/reduction_indices?
"batch_normalization_3/moments/meanMeanconv3d_3/Elu:activations:0=batch_normalization_3/moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2$
"batch_normalization_3/moments/mean?
*batch_normalization_3/moments/StopGradientStopGradient+batch_normalization_3/moments/mean:output:0*
T0**
_output_shapes
:2,
*batch_normalization_3/moments/StopGradient?
/batch_normalization_3/moments/SquaredDifferenceSquaredDifferenceconv3d_3/Elu:activations:03batch_normalization_3/moments/StopGradient:output:0*
T0*3
_output_shapes!
:?????????21
/batch_normalization_3/moments/SquaredDifference?
8batch_normalization_3/moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2:
8batch_normalization_3/moments/variance/reduction_indices?
&batch_normalization_3/moments/varianceMean3batch_normalization_3/moments/SquaredDifference:z:0Abatch_normalization_3/moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2(
&batch_normalization_3/moments/variance?
%batch_normalization_3/moments/SqueezeSqueeze+batch_normalization_3/moments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2'
%batch_normalization_3/moments/Squeeze?
'batch_normalization_3/moments/Squeeze_1Squeeze/batch_normalization_3/moments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2)
'batch_normalization_3/moments/Squeeze_1?
+batch_normalization_3/AssignMovingAvg/decayConst*@
_class6
42loc:@batch_normalization_3/AssignMovingAvg/1226561*
_output_shapes
: *
dtype0*
valueB
 *
?#<2-
+batch_normalization_3/AssignMovingAvg/decay?
4batch_normalization_3/AssignMovingAvg/ReadVariableOpReadVariableOp-batch_normalization_3_assignmovingavg_1226561*
_output_shapes
:*
dtype026
4batch_normalization_3/AssignMovingAvg/ReadVariableOp?
)batch_normalization_3/AssignMovingAvg/subSub<batch_normalization_3/AssignMovingAvg/ReadVariableOp:value:0.batch_normalization_3/moments/Squeeze:output:0*
T0*@
_class6
42loc:@batch_normalization_3/AssignMovingAvg/1226561*
_output_shapes
:2+
)batch_normalization_3/AssignMovingAvg/sub?
)batch_normalization_3/AssignMovingAvg/mulMul-batch_normalization_3/AssignMovingAvg/sub:z:04batch_normalization_3/AssignMovingAvg/decay:output:0*
T0*@
_class6
42loc:@batch_normalization_3/AssignMovingAvg/1226561*
_output_shapes
:2+
)batch_normalization_3/AssignMovingAvg/mul?
9batch_normalization_3/AssignMovingAvg/AssignSubVariableOpAssignSubVariableOp-batch_normalization_3_assignmovingavg_1226561-batch_normalization_3/AssignMovingAvg/mul:z:05^batch_normalization_3/AssignMovingAvg/ReadVariableOp*@
_class6
42loc:@batch_normalization_3/AssignMovingAvg/1226561*
_output_shapes
 *
dtype02;
9batch_normalization_3/AssignMovingAvg/AssignSubVariableOp?
-batch_normalization_3/AssignMovingAvg_1/decayConst*B
_class8
64loc:@batch_normalization_3/AssignMovingAvg_1/1226567*
_output_shapes
: *
dtype0*
valueB
 *
?#<2/
-batch_normalization_3/AssignMovingAvg_1/decay?
6batch_normalization_3/AssignMovingAvg_1/ReadVariableOpReadVariableOp/batch_normalization_3_assignmovingavg_1_1226567*
_output_shapes
:*
dtype028
6batch_normalization_3/AssignMovingAvg_1/ReadVariableOp?
+batch_normalization_3/AssignMovingAvg_1/subSub>batch_normalization_3/AssignMovingAvg_1/ReadVariableOp:value:00batch_normalization_3/moments/Squeeze_1:output:0*
T0*B
_class8
64loc:@batch_normalization_3/AssignMovingAvg_1/1226567*
_output_shapes
:2-
+batch_normalization_3/AssignMovingAvg_1/sub?
+batch_normalization_3/AssignMovingAvg_1/mulMul/batch_normalization_3/AssignMovingAvg_1/sub:z:06batch_normalization_3/AssignMovingAvg_1/decay:output:0*
T0*B
_class8
64loc:@batch_normalization_3/AssignMovingAvg_1/1226567*
_output_shapes
:2-
+batch_normalization_3/AssignMovingAvg_1/mul?
;batch_normalization_3/AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOp/batch_normalization_3_assignmovingavg_1_1226567/batch_normalization_3/AssignMovingAvg_1/mul:z:07^batch_normalization_3/AssignMovingAvg_1/ReadVariableOp*B
_class8
64loc:@batch_normalization_3/AssignMovingAvg_1/1226567*
_output_shapes
 *
dtype02=
;batch_normalization_3/AssignMovingAvg_1/AssignSubVariableOp?
%batch_normalization_3/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2'
%batch_normalization_3/batchnorm/add/y?
#batch_normalization_3/batchnorm/addAddV20batch_normalization_3/moments/Squeeze_1:output:0.batch_normalization_3/batchnorm/add/y:output:0*
T0*
_output_shapes
:2%
#batch_normalization_3/batchnorm/add?
%batch_normalization_3/batchnorm/RsqrtRsqrt'batch_normalization_3/batchnorm/add:z:0*
T0*
_output_shapes
:2'
%batch_normalization_3/batchnorm/Rsqrt?
2batch_normalization_3/batchnorm/mul/ReadVariableOpReadVariableOp;batch_normalization_3_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype024
2batch_normalization_3/batchnorm/mul/ReadVariableOp?
#batch_normalization_3/batchnorm/mulMul)batch_normalization_3/batchnorm/Rsqrt:y:0:batch_normalization_3/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2%
#batch_normalization_3/batchnorm/mul?
%batch_normalization_3/batchnorm/mul_1Mulconv3d_3/Elu:activations:0'batch_normalization_3/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2'
%batch_normalization_3/batchnorm/mul_1?
%batch_normalization_3/batchnorm/mul_2Mul.batch_normalization_3/moments/Squeeze:output:0'batch_normalization_3/batchnorm/mul:z:0*
T0*
_output_shapes
:2'
%batch_normalization_3/batchnorm/mul_2?
.batch_normalization_3/batchnorm/ReadVariableOpReadVariableOp7batch_normalization_3_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype020
.batch_normalization_3/batchnorm/ReadVariableOp?
#batch_normalization_3/batchnorm/subSub6batch_normalization_3/batchnorm/ReadVariableOp:value:0)batch_normalization_3/batchnorm/mul_2:z:0*
T0*
_output_shapes
:2%
#batch_normalization_3/batchnorm/sub?
%batch_normalization_3/batchnorm/add_1AddV2)batch_normalization_3/batchnorm/mul_1:z:0'batch_normalization_3/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2'
%batch_normalization_3/batchnorm/add_1?
average_pooling3d/AvgPool3D	AvgPool3D)batch_normalization_3/batchnorm/add_1:z:0*
T0*3
_output_shapes!
:?????????*
ksize	
*
paddingVALID*
strides	
2
average_pooling3d/AvgPool3Do
flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"????   2
flatten/Const?
flatten/ReshapeReshape$average_pooling3d/AvgPool3D:output:0flatten/Const:output:0*
T0*(
_output_shapes
:??????????
2
flatten/Reshapes
dropout/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *UU??2
dropout/dropout/Const?
dropout/dropout/MulMulflatten/Reshape:output:0dropout/dropout/Const:output:0*
T0*(
_output_shapes
:??????????
2
dropout/dropout/Mulv
dropout/dropout/ShapeShapeflatten/Reshape:output:0*
T0*
_output_shapes
:2
dropout/dropout/Shape?
,dropout/dropout/random_uniform/RandomUniformRandomUniformdropout/dropout/Shape:output:0*
T0*(
_output_shapes
:??????????
*
dtype02.
,dropout/dropout/random_uniform/RandomUniform?
dropout/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2 
dropout/dropout/GreaterEqual/y?
dropout/dropout/GreaterEqualGreaterEqual5dropout/dropout/random_uniform/RandomUniform:output:0'dropout/dropout/GreaterEqual/y:output:0*
T0*(
_output_shapes
:??????????
2
dropout/dropout/GreaterEqual?
dropout/dropout/CastCast dropout/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*(
_output_shapes
:??????????
2
dropout/dropout/Cast?
dropout/dropout/Mul_1Muldropout/dropout/Mul:z:0dropout/dropout/Cast:y:0*
T0*(
_output_shapes
:??????????
2
dropout/dropout/Mul_1?
layer1/MatMul/ReadVariableOpReadVariableOp%layer1_matmul_readvariableop_resource* 
_output_shapes
:
?
?*
dtype02
layer1/MatMul/ReadVariableOp?
layer1/MatMulMatMuldropout/dropout/Mul_1:z:0$layer1/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2
layer1/MatMul?
layer1/BiasAdd/ReadVariableOpReadVariableOp&layer1_biasadd_readvariableop_resource*
_output_shapes	
:?*
dtype02
layer1/BiasAdd/ReadVariableOp?
layer1/BiasAddBiasAddlayer1/MatMul:product:0%layer1/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2
layer1/BiasAddk

layer1/EluElulayer1/BiasAdd:output:0*
T0*(
_output_shapes
:??????????2

layer1/Elu?
IdentityIdentitylayer1/Elu:activations:08^batch_normalization/AssignMovingAvg/AssignSubVariableOp:^batch_normalization/AssignMovingAvg_1/AssignSubVariableOp:^batch_normalization_1/AssignMovingAvg/AssignSubVariableOp<^batch_normalization_1/AssignMovingAvg_1/AssignSubVariableOp:^batch_normalization_2/AssignMovingAvg/AssignSubVariableOp<^batch_normalization_2/AssignMovingAvg_1/AssignSubVariableOp:^batch_normalization_3/AssignMovingAvg/AssignSubVariableOp<^batch_normalization_3/AssignMovingAvg_1/AssignSubVariableOp*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????::::::::::::::::::::::::::2r
7batch_normalization/AssignMovingAvg/AssignSubVariableOp7batch_normalization/AssignMovingAvg/AssignSubVariableOp2v
9batch_normalization/AssignMovingAvg_1/AssignSubVariableOp9batch_normalization/AssignMovingAvg_1/AssignSubVariableOp2v
9batch_normalization_1/AssignMovingAvg/AssignSubVariableOp9batch_normalization_1/AssignMovingAvg/AssignSubVariableOp2z
;batch_normalization_1/AssignMovingAvg_1/AssignSubVariableOp;batch_normalization_1/AssignMovingAvg_1/AssignSubVariableOp2v
9batch_normalization_2/AssignMovingAvg/AssignSubVariableOp9batch_normalization_2/AssignMovingAvg/AssignSubVariableOp2z
;batch_normalization_2/AssignMovingAvg_1/AssignSubVariableOp;batch_normalization_2/AssignMovingAvg_1/AssignSubVariableOp2v
9batch_normalization_3/AssignMovingAvg/AssignSubVariableOp9batch_normalization_3/AssignMovingAvg/AssignSubVariableOp2z
;batch_normalization_3/AssignMovingAvg_1/AssignSubVariableOp;batch_normalization_3/AssignMovingAvg_1/AssignSubVariableOp:\ X
4
_output_shapes"
 :??????????
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_1222601

inputs%
!batchnorm_readvariableop_resource)
%batchnorm_mul_readvariableop_resource'
#batchnorm_readvariableop_1_resource'
#batchnorm_readvariableop_2_resource
identity??
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/mul_1?
batchnorm/ReadVariableOp_1ReadVariableOp#batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_1?
batchnorm/mul_2Mul"batchnorm/ReadVariableOp_1:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOp_2ReadVariableOp#batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_2?
batchnorm/subSub"batchnorm/ReadVariableOp_2:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????:::::v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_1222741

inputs%
!batchnorm_readvariableop_resource)
%batchnorm_mul_readvariableop_resource'
#batchnorm_readvariableop_1_resource'
#batchnorm_readvariableop_2_resource
identity??
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/mul_1?
batchnorm/ReadVariableOp_1ReadVariableOp#batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_1?
batchnorm/mul_2Mul"batchnorm/ReadVariableOp_1:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOp_2ReadVariableOp#batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_2?
batchnorm/subSub"batchnorm/ReadVariableOp_2:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????:::::v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
?
?
,__inference_sequential_layer_call_fn_1224105
functional_1_input
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16

unknown_17

unknown_18

unknown_19

unknown_20

unknown_21

unknown_22

unknown_23

unknown_24

unknown_25

unknown_26

unknown_27

unknown_28
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallfunctional_1_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19
unknown_20
unknown_21
unknown_22
unknown_23
unknown_24
unknown_25
unknown_26
unknown_27
unknown_28**
Tin#
!2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
*8
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *P
fKRI
G__inference_sequential_layer_call_and_return_conditional_losses_12240422
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????
2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????::::::::::::::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:h d
4
_output_shapes"
 :??????????
,
_user_specified_namefunctional_1_input
?
?
5__inference_batch_normalization_layer_call_fn_1227080

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *Y
fTRR
P__inference_batch_normalization_layer_call_and_return_conditional_losses_12228292
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????::::22
StatefulPartitionedCallStatefulPartitionedCall:[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?
?
7__inference_batch_normalization_1_layer_call_fn_1227264

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_12229472
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????::::22
StatefulPartitionedCallStatefulPartitionedCall:[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?
?
P__inference_batch_normalization_layer_call_and_return_conditional_losses_1226985

inputs%
!batchnorm_readvariableop_resource)
%batchnorm_mul_readvariableop_resource'
#batchnorm_readvariableop_1_resource'
#batchnorm_readvariableop_2_resource
identity??
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/mul_1?
batchnorm/ReadVariableOp_1ReadVariableOp#batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_1?
batchnorm/mul_2Mul"batchnorm/ReadVariableOp_1:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOp_2ReadVariableOp#batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_2?
batchnorm/subSub"batchnorm/ReadVariableOp_2:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????:::::v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
?C
?	
I__inference_functional_1_layer_call_and_return_conditional_losses_1223312
input_1
conv3d_1222789
conv3d_1222791
batch_normalization_1222876
batch_normalization_1222878
batch_normalization_1222880
batch_normalization_1222882
conv3d_1_1222907
conv3d_1_1222909!
batch_normalization_1_1222994!
batch_normalization_1_1222996!
batch_normalization_1_1222998!
batch_normalization_1_1223000
conv3d_2_1223025
conv3d_2_1223027!
batch_normalization_2_1223112!
batch_normalization_2_1223114!
batch_normalization_2_1223116!
batch_normalization_2_1223118
conv3d_3_1223143
conv3d_3_1223145!
batch_normalization_3_1223230!
batch_normalization_3_1223232!
batch_normalization_3_1223234!
batch_normalization_3_1223236
layer1_1223306
layer1_1223308
identity??+batch_normalization/StatefulPartitionedCall?-batch_normalization_1/StatefulPartitionedCall?-batch_normalization_2/StatefulPartitionedCall?-batch_normalization_3/StatefulPartitionedCall?conv3d/StatefulPartitionedCall? conv3d_1/StatefulPartitionedCall? conv3d_2/StatefulPartitionedCall? conv3d_3/StatefulPartitionedCall?dropout/StatefulPartitionedCall?layer1/StatefulPartitionedCall?
conv3d/StatefulPartitionedCallStatefulPartitionedCallinput_1conv3d_1222789conv3d_1222791*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *L
fGRE
C__inference_conv3d_layer_call_and_return_conditional_losses_12227782 
conv3d/StatefulPartitionedCall?
+batch_normalization/StatefulPartitionedCallStatefulPartitionedCall'conv3d/StatefulPartitionedCall:output:0batch_normalization_1222876batch_normalization_1222878batch_normalization_1222880batch_normalization_1222882*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *Y
fTRR
P__inference_batch_normalization_layer_call_and_return_conditional_losses_12228292-
+batch_normalization/StatefulPartitionedCall?
 conv3d_1/StatefulPartitionedCallStatefulPartitionedCall4batch_normalization/StatefulPartitionedCall:output:0conv3d_1_1222907conv3d_1_1222909*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv3d_1_layer_call_and_return_conditional_losses_12228962"
 conv3d_1/StatefulPartitionedCall?
-batch_normalization_1/StatefulPartitionedCallStatefulPartitionedCall)conv3d_1/StatefulPartitionedCall:output:0batch_normalization_1_1222994batch_normalization_1_1222996batch_normalization_1_1222998batch_normalization_1_1223000*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_12229472/
-batch_normalization_1/StatefulPartitionedCall?
 conv3d_2/StatefulPartitionedCallStatefulPartitionedCall6batch_normalization_1/StatefulPartitionedCall:output:0conv3d_2_1223025conv3d_2_1223027*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv3d_2_layer_call_and_return_conditional_losses_12230142"
 conv3d_2/StatefulPartitionedCall?
-batch_normalization_2/StatefulPartitionedCallStatefulPartitionedCall)conv3d_2/StatefulPartitionedCall:output:0batch_normalization_2_1223112batch_normalization_2_1223114batch_normalization_2_1223116batch_normalization_2_1223118*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_12230652/
-batch_normalization_2/StatefulPartitionedCall?
 conv3d_3/StatefulPartitionedCallStatefulPartitionedCall6batch_normalization_2/StatefulPartitionedCall:output:0conv3d_3_1223143conv3d_3_1223145*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv3d_3_layer_call_and_return_conditional_losses_12231322"
 conv3d_3/StatefulPartitionedCall?
-batch_normalization_3/StatefulPartitionedCallStatefulPartitionedCall)conv3d_3/StatefulPartitionedCall:output:0batch_normalization_3_1223230batch_normalization_3_1223232batch_normalization_3_1223234batch_normalization_3_1223236*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_12231832/
-batch_normalization_3/StatefulPartitionedCall?
!average_pooling3d/PartitionedCallPartitionedCall6batch_normalization_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *W
fRRP
N__inference_average_pooling3d_layer_call_and_return_conditional_losses_12227582#
!average_pooling3d/PartitionedCall?
flatten/PartitionedCallPartitionedCall*average_pooling3d/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????
* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_flatten_layer_call_and_return_conditional_losses_12232462
flatten/PartitionedCall?
dropout/StatefulPartitionedCallStatefulPartitionedCall flatten/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????
* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dropout_layer_call_and_return_conditional_losses_12232662!
dropout/StatefulPartitionedCall?
layer1/StatefulPartitionedCallStatefulPartitionedCall(dropout/StatefulPartitionedCall:output:0layer1_1223306layer1_1223308*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *L
fGRE
C__inference_layer1_layer_call_and_return_conditional_losses_12232952 
layer1/StatefulPartitionedCall?
IdentityIdentity'layer1/StatefulPartitionedCall:output:0,^batch_normalization/StatefulPartitionedCall.^batch_normalization_1/StatefulPartitionedCall.^batch_normalization_2/StatefulPartitionedCall.^batch_normalization_3/StatefulPartitionedCall^conv3d/StatefulPartitionedCall!^conv3d_1/StatefulPartitionedCall!^conv3d_2/StatefulPartitionedCall!^conv3d_3/StatefulPartitionedCall ^dropout/StatefulPartitionedCall^layer1/StatefulPartitionedCall*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????::::::::::::::::::::::::::2Z
+batch_normalization/StatefulPartitionedCall+batch_normalization/StatefulPartitionedCall2^
-batch_normalization_1/StatefulPartitionedCall-batch_normalization_1/StatefulPartitionedCall2^
-batch_normalization_2/StatefulPartitionedCall-batch_normalization_2/StatefulPartitionedCall2^
-batch_normalization_3/StatefulPartitionedCall-batch_normalization_3/StatefulPartitionedCall2@
conv3d/StatefulPartitionedCallconv3d/StatefulPartitionedCall2D
 conv3d_1/StatefulPartitionedCall conv3d_1/StatefulPartitionedCall2D
 conv3d_2/StatefulPartitionedCall conv3d_2/StatefulPartitionedCall2D
 conv3d_3/StatefulPartitionedCall conv3d_3/StatefulPartitionedCall2B
dropout/StatefulPartitionedCalldropout/StatefulPartitionedCall2@
layer1/StatefulPartitionedCalllayer1/StatefulPartitionedCall:] Y
4
_output_shapes"
 :??????????
!
_user_specified_name	input_1
?
E
)__inference_flatten_layer_call_fn_1227656

inputs
identity?
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????
* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_flatten_layer_call_and_return_conditional_losses_12232462
PartitionedCallm
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:??????????
2

Identity"
identityIdentity:output:0*2
_input_shapes!
:?????????:[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?*
?
P__inference_batch_normalization_layer_call_and_return_conditional_losses_1227047

inputs
assignmovingavg_1227022
assignmovingavg_1_1227028)
%batchnorm_mul_readvariableop_resource%
!batchnorm_readvariableop_resource
identity??#AssignMovingAvg/AssignSubVariableOp?%AssignMovingAvg_1/AssignSubVariableOp?
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2 
moments/mean/reduction_indices?
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/mean?
moments/StopGradientStopGradientmoments/mean:output:0*
T0**
_output_shapes
:2
moments/StopGradient?
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*3
_output_shapes!
:?????????2
moments/SquaredDifference?
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2$
"moments/variance/reduction_indices?
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/variance?
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze?
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze_1?
AssignMovingAvg/decayConst**
_class 
loc:@AssignMovingAvg/1227022*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg/decay?
AssignMovingAvg/ReadVariableOpReadVariableOpassignmovingavg_1227022*
_output_shapes
:*
dtype02 
AssignMovingAvg/ReadVariableOp?
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0**
_class 
loc:@AssignMovingAvg/1227022*
_output_shapes
:2
AssignMovingAvg/sub?
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0**
_class 
loc:@AssignMovingAvg/1227022*
_output_shapes
:2
AssignMovingAvg/mul?
#AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1227022AssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp**
_class 
loc:@AssignMovingAvg/1227022*
_output_shapes
 *
dtype02%
#AssignMovingAvg/AssignSubVariableOp?
AssignMovingAvg_1/decayConst*,
_class"
 loc:@AssignMovingAvg_1/1227028*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg_1/decay?
 AssignMovingAvg_1/ReadVariableOpReadVariableOpassignmovingavg_1_1227028*
_output_shapes
:*
dtype02"
 AssignMovingAvg_1/ReadVariableOp?
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1227028*
_output_shapes
:2
AssignMovingAvg_1/sub?
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1227028*
_output_shapes
:2
AssignMovingAvg_1/mul?
%AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1_1227028AssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*,
_class"
 loc:@AssignMovingAvg_1/1227028*
_output_shapes
 *
dtype02'
%AssignMovingAvg_1/AssignSubVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/mul_1{
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp?
batchnorm/subSub batchnorm/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0$^AssignMovingAvg/AssignSubVariableOp&^AssignMovingAvg_1/AssignSubVariableOp*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????::::2J
#AssignMovingAvg/AssignSubVariableOp#AssignMovingAvg/AssignSubVariableOp2N
%AssignMovingAvg_1/AssignSubVariableOp%AssignMovingAvg_1/AssignSubVariableOp:[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?

o
__inference_loss_fn_1_1226910=
9dense_1_kernel_regularizer_square_readvariableop_resource
identity??
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp9dense_1_kernel_regularizer_square_readvariableop_resource*
_output_shapes

:d
*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp?
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d
2#
!dense_1/kernel/Regularizer/Square?
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const?
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/Sum?
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2"
 dense_1/kernel/Regularizer/mul/x?
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/mule
IdentityIdentity"dense_1/kernel/Regularizer/mul:z:0*
T0*
_output_shapes
: 2

Identity"
identityIdentity:output:0*
_input_shapes
:
?,
?
G__inference_sequential_layer_call_and_return_conditional_losses_1223881
functional_1_input
functional_1_1223749
functional_1_1223751
functional_1_1223753
functional_1_1223755
functional_1_1223757
functional_1_1223759
functional_1_1223761
functional_1_1223763
functional_1_1223765
functional_1_1223767
functional_1_1223769
functional_1_1223771
functional_1_1223773
functional_1_1223775
functional_1_1223777
functional_1_1223779
functional_1_1223781
functional_1_1223783
functional_1_1223785
functional_1_1223787
functional_1_1223789
functional_1_1223791
functional_1_1223793
functional_1_1223795
functional_1_1223797
functional_1_1223799
dense_1223830
dense_1223832
dense_1_1223863
dense_1_1223865
identity??dense/StatefulPartitionedCall?dense_1/StatefulPartitionedCall?$functional_1/StatefulPartitionedCall?
$functional_1/StatefulPartitionedCallStatefulPartitionedCallfunctional_1_inputfunctional_1_1223749functional_1_1223751functional_1_1223753functional_1_1223755functional_1_1223757functional_1_1223759functional_1_1223761functional_1_1223763functional_1_1223765functional_1_1223767functional_1_1223769functional_1_1223771functional_1_1223773functional_1_1223775functional_1_1223777functional_1_1223779functional_1_1223781functional_1_1223783functional_1_1223785functional_1_1223787functional_1_1223789functional_1_1223791functional_1_1223793functional_1_1223795functional_1_1223797functional_1_1223799*&
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*4
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *R
fMRK
I__inference_functional_1_layer_call_and_return_conditional_losses_12234512&
$functional_1/StatefulPartitionedCall?
dense/StatefulPartitionedCallStatefulPartitionedCall-functional_1/StatefulPartitionedCall:output:0dense_1223830dense_1223832*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????d*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *K
fFRD
B__inference_dense_layer_call_and_return_conditional_losses_12238192
dense/StatefulPartitionedCall?
dense_1/StatefulPartitionedCallStatefulPartitionedCall&dense/StatefulPartitionedCall:output:0dense_1_1223863dense_1_1223865*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_1_layer_call_and_return_conditional_losses_12238522!
dense_1/StatefulPartitionedCall?
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_1223830*
_output_shapes
:	?d*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOp?
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?d2!
dense/kernel/Regularizer/Square?
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const?
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/Sum?
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2 
dense/kernel/Regularizer/mul/x?
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/mul?
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_1_1223863*
_output_shapes

:d
*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp?
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d
2#
!dense_1/kernel/Regularizer/Square?
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const?
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/Sum?
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2"
 dense_1/kernel/Regularizer/mul/x?
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/mul?
IdentityIdentity(dense_1/StatefulPartitionedCall:output:0^dense/StatefulPartitionedCall ^dense_1/StatefulPartitionedCall%^functional_1/StatefulPartitionedCall*
T0*'
_output_shapes
:?????????
2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????::::::::::::::::::::::::::::::2>
dense/StatefulPartitionedCalldense/StatefulPartitionedCall2B
dense_1/StatefulPartitionedCalldense_1/StatefulPartitionedCall2L
$functional_1/StatefulPartitionedCall$functional_1/StatefulPartitionedCall:h d
4
_output_shapes"
 :??????????
,
_user_specified_namefunctional_1_input
?*
?
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_1227599

inputs
assignmovingavg_1227574
assignmovingavg_1_1227580)
%batchnorm_mul_readvariableop_resource%
!batchnorm_readvariableop_resource
identity??#AssignMovingAvg/AssignSubVariableOp?%AssignMovingAvg_1/AssignSubVariableOp?
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2 
moments/mean/reduction_indices?
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/mean?
moments/StopGradientStopGradientmoments/mean:output:0*
T0**
_output_shapes
:2
moments/StopGradient?
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*3
_output_shapes!
:?????????2
moments/SquaredDifference?
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2$
"moments/variance/reduction_indices?
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/variance?
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze?
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze_1?
AssignMovingAvg/decayConst**
_class 
loc:@AssignMovingAvg/1227574*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg/decay?
AssignMovingAvg/ReadVariableOpReadVariableOpassignmovingavg_1227574*
_output_shapes
:*
dtype02 
AssignMovingAvg/ReadVariableOp?
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0**
_class 
loc:@AssignMovingAvg/1227574*
_output_shapes
:2
AssignMovingAvg/sub?
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0**
_class 
loc:@AssignMovingAvg/1227574*
_output_shapes
:2
AssignMovingAvg/mul?
#AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1227574AssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp**
_class 
loc:@AssignMovingAvg/1227574*
_output_shapes
 *
dtype02%
#AssignMovingAvg/AssignSubVariableOp?
AssignMovingAvg_1/decayConst*,
_class"
 loc:@AssignMovingAvg_1/1227580*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg_1/decay?
 AssignMovingAvg_1/ReadVariableOpReadVariableOpassignmovingavg_1_1227580*
_output_shapes
:*
dtype02"
 AssignMovingAvg_1/ReadVariableOp?
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1227580*
_output_shapes
:2
AssignMovingAvg_1/sub?
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1227580*
_output_shapes
:2
AssignMovingAvg_1/mul?
%AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1_1227580AssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*,
_class"
 loc:@AssignMovingAvg_1/1227580*
_output_shapes
 *
dtype02'
%AssignMovingAvg_1/AssignSubVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/mul_1{
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp?
batchnorm/subSub batchnorm/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0$^AssignMovingAvg/AssignSubVariableOp&^AssignMovingAvg_1/AssignSubVariableOp*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????::::2J
#AssignMovingAvg/AssignSubVariableOp#AssignMovingAvg/AssignSubVariableOp2N
%AssignMovingAvg_1/AssignSubVariableOp%AssignMovingAvg_1/AssignSubVariableOp:[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?
?
7__inference_batch_normalization_3_layer_call_fn_1227550

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *N
_output_shapes<
::8????????????????????????????????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_12227082
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????::::22
StatefulPartitionedCallStatefulPartitionedCall:v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
?
?
P__inference_batch_normalization_layer_call_and_return_conditional_losses_1222321

inputs%
!batchnorm_readvariableop_resource)
%batchnorm_mul_readvariableop_resource'
#batchnorm_readvariableop_1_resource'
#batchnorm_readvariableop_2_resource
identity??
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/mul_1?
batchnorm/ReadVariableOp_1ReadVariableOp#batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_1?
batchnorm/mul_2Mul"batchnorm/ReadVariableOp_1:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOp_2ReadVariableOp#batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_2?
batchnorm/subSub"batchnorm/ReadVariableOp_2:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????:::::v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_1227251

inputs%
!batchnorm_readvariableop_resource)
%batchnorm_mul_readvariableop_resource'
#batchnorm_readvariableop_1_resource'
#batchnorm_readvariableop_2_resource
identity??
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/mul_1?
batchnorm/ReadVariableOp_1ReadVariableOp#batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_1?
batchnorm/mul_2Mul"batchnorm/ReadVariableOp_1:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOp_2ReadVariableOp#batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_2?
batchnorm/subSub"batchnorm/ReadVariableOp_2:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/add_1s
IdentityIdentitybatchnorm/add_1:z:0*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????:::::[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?	
?
E__inference_conv3d_1_layer_call_and_return_conditional_losses_1222896

inputs"
conv3d_readvariableop_resource#
biasadd_readvariableop_resource
identity??
Conv3D/ReadVariableOpReadVariableOpconv3d_readvariableop_resource**
_output_shapes
:*
dtype02
Conv3D/ReadVariableOp?
Conv3DConv3DinputsConv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
Conv3D?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddConv3D:output:0BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2	
BiasAdda
EluEluBiasAdd:output:0*
T0*3
_output_shapes!
:?????????2
Eluq
IdentityIdentityElu:activations:0*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*:
_input_shapes)
':?????????:::[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
??
?
G__inference_sequential_layer_call_and_return_conditional_losses_1226107

inputs6
2functional_1_conv3d_conv3d_readvariableop_resource7
3functional_1_conv3d_biasadd_readvariableop_resource<
8functional_1_batch_normalization_assignmovingavg_1225921>
:functional_1_batch_normalization_assignmovingavg_1_1225927J
Ffunctional_1_batch_normalization_batchnorm_mul_readvariableop_resourceF
Bfunctional_1_batch_normalization_batchnorm_readvariableop_resource8
4functional_1_conv3d_1_conv3d_readvariableop_resource9
5functional_1_conv3d_1_biasadd_readvariableop_resource>
:functional_1_batch_normalization_1_assignmovingavg_1225960@
<functional_1_batch_normalization_1_assignmovingavg_1_1225966L
Hfunctional_1_batch_normalization_1_batchnorm_mul_readvariableop_resourceH
Dfunctional_1_batch_normalization_1_batchnorm_readvariableop_resource8
4functional_1_conv3d_2_conv3d_readvariableop_resource9
5functional_1_conv3d_2_biasadd_readvariableop_resource>
:functional_1_batch_normalization_2_assignmovingavg_1225999@
<functional_1_batch_normalization_2_assignmovingavg_1_1226005L
Hfunctional_1_batch_normalization_2_batchnorm_mul_readvariableop_resourceH
Dfunctional_1_batch_normalization_2_batchnorm_readvariableop_resource8
4functional_1_conv3d_3_conv3d_readvariableop_resource9
5functional_1_conv3d_3_biasadd_readvariableop_resource>
:functional_1_batch_normalization_3_assignmovingavg_1226038@
<functional_1_batch_normalization_3_assignmovingavg_1_1226044L
Hfunctional_1_batch_normalization_3_batchnorm_mul_readvariableop_resourceH
Dfunctional_1_batch_normalization_3_batchnorm_readvariableop_resource6
2functional_1_layer1_matmul_readvariableop_resource7
3functional_1_layer1_biasadd_readvariableop_resource(
$dense_matmul_readvariableop_resource)
%dense_biasadd_readvariableop_resource*
&dense_1_matmul_readvariableop_resource+
'dense_1_biasadd_readvariableop_resource
identity??Dfunctional_1/batch_normalization/AssignMovingAvg/AssignSubVariableOp?Ffunctional_1/batch_normalization/AssignMovingAvg_1/AssignSubVariableOp?Ffunctional_1/batch_normalization_1/AssignMovingAvg/AssignSubVariableOp?Hfunctional_1/batch_normalization_1/AssignMovingAvg_1/AssignSubVariableOp?Ffunctional_1/batch_normalization_2/AssignMovingAvg/AssignSubVariableOp?Hfunctional_1/batch_normalization_2/AssignMovingAvg_1/AssignSubVariableOp?Ffunctional_1/batch_normalization_3/AssignMovingAvg/AssignSubVariableOp?Hfunctional_1/batch_normalization_3/AssignMovingAvg_1/AssignSubVariableOp?
)functional_1/conv3d/Conv3D/ReadVariableOpReadVariableOp2functional_1_conv3d_conv3d_readvariableop_resource*+
_output_shapes
:?*
dtype02+
)functional_1/conv3d/Conv3D/ReadVariableOp?
functional_1/conv3d/Conv3DConv3Dinputs1functional_1/conv3d/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
functional_1/conv3d/Conv3D?
*functional_1/conv3d/BiasAdd/ReadVariableOpReadVariableOp3functional_1_conv3d_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02,
*functional_1/conv3d/BiasAdd/ReadVariableOp?
functional_1/conv3d/BiasAddBiasAdd#functional_1/conv3d/Conv3D:output:02functional_1/conv3d/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2
functional_1/conv3d/BiasAdd?
?functional_1/batch_normalization/moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2A
?functional_1/batch_normalization/moments/mean/reduction_indices?
-functional_1/batch_normalization/moments/meanMean$functional_1/conv3d/BiasAdd:output:0Hfunctional_1/batch_normalization/moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2/
-functional_1/batch_normalization/moments/mean?
5functional_1/batch_normalization/moments/StopGradientStopGradient6functional_1/batch_normalization/moments/mean:output:0*
T0**
_output_shapes
:27
5functional_1/batch_normalization/moments/StopGradient?
:functional_1/batch_normalization/moments/SquaredDifferenceSquaredDifference$functional_1/conv3d/BiasAdd:output:0>functional_1/batch_normalization/moments/StopGradient:output:0*
T0*3
_output_shapes!
:?????????2<
:functional_1/batch_normalization/moments/SquaredDifference?
Cfunctional_1/batch_normalization/moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2E
Cfunctional_1/batch_normalization/moments/variance/reduction_indices?
1functional_1/batch_normalization/moments/varianceMean>functional_1/batch_normalization/moments/SquaredDifference:z:0Lfunctional_1/batch_normalization/moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(23
1functional_1/batch_normalization/moments/variance?
0functional_1/batch_normalization/moments/SqueezeSqueeze6functional_1/batch_normalization/moments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 22
0functional_1/batch_normalization/moments/Squeeze?
2functional_1/batch_normalization/moments/Squeeze_1Squeeze:functional_1/batch_normalization/moments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 24
2functional_1/batch_normalization/moments/Squeeze_1?
6functional_1/batch_normalization/AssignMovingAvg/decayConst*K
_classA
?=loc:@functional_1/batch_normalization/AssignMovingAvg/1225921*
_output_shapes
: *
dtype0*
valueB
 *
?#<28
6functional_1/batch_normalization/AssignMovingAvg/decay?
?functional_1/batch_normalization/AssignMovingAvg/ReadVariableOpReadVariableOp8functional_1_batch_normalization_assignmovingavg_1225921*
_output_shapes
:*
dtype02A
?functional_1/batch_normalization/AssignMovingAvg/ReadVariableOp?
4functional_1/batch_normalization/AssignMovingAvg/subSubGfunctional_1/batch_normalization/AssignMovingAvg/ReadVariableOp:value:09functional_1/batch_normalization/moments/Squeeze:output:0*
T0*K
_classA
?=loc:@functional_1/batch_normalization/AssignMovingAvg/1225921*
_output_shapes
:26
4functional_1/batch_normalization/AssignMovingAvg/sub?
4functional_1/batch_normalization/AssignMovingAvg/mulMul8functional_1/batch_normalization/AssignMovingAvg/sub:z:0?functional_1/batch_normalization/AssignMovingAvg/decay:output:0*
T0*K
_classA
?=loc:@functional_1/batch_normalization/AssignMovingAvg/1225921*
_output_shapes
:26
4functional_1/batch_normalization/AssignMovingAvg/mul?
Dfunctional_1/batch_normalization/AssignMovingAvg/AssignSubVariableOpAssignSubVariableOp8functional_1_batch_normalization_assignmovingavg_12259218functional_1/batch_normalization/AssignMovingAvg/mul:z:0@^functional_1/batch_normalization/AssignMovingAvg/ReadVariableOp*K
_classA
?=loc:@functional_1/batch_normalization/AssignMovingAvg/1225921*
_output_shapes
 *
dtype02F
Dfunctional_1/batch_normalization/AssignMovingAvg/AssignSubVariableOp?
8functional_1/batch_normalization/AssignMovingAvg_1/decayConst*M
_classC
A?loc:@functional_1/batch_normalization/AssignMovingAvg_1/1225927*
_output_shapes
: *
dtype0*
valueB
 *
?#<2:
8functional_1/batch_normalization/AssignMovingAvg_1/decay?
Afunctional_1/batch_normalization/AssignMovingAvg_1/ReadVariableOpReadVariableOp:functional_1_batch_normalization_assignmovingavg_1_1225927*
_output_shapes
:*
dtype02C
Afunctional_1/batch_normalization/AssignMovingAvg_1/ReadVariableOp?
6functional_1/batch_normalization/AssignMovingAvg_1/subSubIfunctional_1/batch_normalization/AssignMovingAvg_1/ReadVariableOp:value:0;functional_1/batch_normalization/moments/Squeeze_1:output:0*
T0*M
_classC
A?loc:@functional_1/batch_normalization/AssignMovingAvg_1/1225927*
_output_shapes
:28
6functional_1/batch_normalization/AssignMovingAvg_1/sub?
6functional_1/batch_normalization/AssignMovingAvg_1/mulMul:functional_1/batch_normalization/AssignMovingAvg_1/sub:z:0Afunctional_1/batch_normalization/AssignMovingAvg_1/decay:output:0*
T0*M
_classC
A?loc:@functional_1/batch_normalization/AssignMovingAvg_1/1225927*
_output_shapes
:28
6functional_1/batch_normalization/AssignMovingAvg_1/mul?
Ffunctional_1/batch_normalization/AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOp:functional_1_batch_normalization_assignmovingavg_1_1225927:functional_1/batch_normalization/AssignMovingAvg_1/mul:z:0B^functional_1/batch_normalization/AssignMovingAvg_1/ReadVariableOp*M
_classC
A?loc:@functional_1/batch_normalization/AssignMovingAvg_1/1225927*
_output_shapes
 *
dtype02H
Ffunctional_1/batch_normalization/AssignMovingAvg_1/AssignSubVariableOp?
0functional_1/batch_normalization/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:22
0functional_1/batch_normalization/batchnorm/add/y?
.functional_1/batch_normalization/batchnorm/addAddV2;functional_1/batch_normalization/moments/Squeeze_1:output:09functional_1/batch_normalization/batchnorm/add/y:output:0*
T0*
_output_shapes
:20
.functional_1/batch_normalization/batchnorm/add?
0functional_1/batch_normalization/batchnorm/RsqrtRsqrt2functional_1/batch_normalization/batchnorm/add:z:0*
T0*
_output_shapes
:22
0functional_1/batch_normalization/batchnorm/Rsqrt?
=functional_1/batch_normalization/batchnorm/mul/ReadVariableOpReadVariableOpFfunctional_1_batch_normalization_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02?
=functional_1/batch_normalization/batchnorm/mul/ReadVariableOp?
.functional_1/batch_normalization/batchnorm/mulMul4functional_1/batch_normalization/batchnorm/Rsqrt:y:0Efunctional_1/batch_normalization/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:20
.functional_1/batch_normalization/batchnorm/mul?
0functional_1/batch_normalization/batchnorm/mul_1Mul$functional_1/conv3d/BiasAdd:output:02functional_1/batch_normalization/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????22
0functional_1/batch_normalization/batchnorm/mul_1?
0functional_1/batch_normalization/batchnorm/mul_2Mul9functional_1/batch_normalization/moments/Squeeze:output:02functional_1/batch_normalization/batchnorm/mul:z:0*
T0*
_output_shapes
:22
0functional_1/batch_normalization/batchnorm/mul_2?
9functional_1/batch_normalization/batchnorm/ReadVariableOpReadVariableOpBfunctional_1_batch_normalization_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02;
9functional_1/batch_normalization/batchnorm/ReadVariableOp?
.functional_1/batch_normalization/batchnorm/subSubAfunctional_1/batch_normalization/batchnorm/ReadVariableOp:value:04functional_1/batch_normalization/batchnorm/mul_2:z:0*
T0*
_output_shapes
:20
.functional_1/batch_normalization/batchnorm/sub?
0functional_1/batch_normalization/batchnorm/add_1AddV24functional_1/batch_normalization/batchnorm/mul_1:z:02functional_1/batch_normalization/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????22
0functional_1/batch_normalization/batchnorm/add_1?
+functional_1/conv3d_1/Conv3D/ReadVariableOpReadVariableOp4functional_1_conv3d_1_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02-
+functional_1/conv3d_1/Conv3D/ReadVariableOp?
functional_1/conv3d_1/Conv3DConv3D4functional_1/batch_normalization/batchnorm/add_1:z:03functional_1/conv3d_1/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
functional_1/conv3d_1/Conv3D?
,functional_1/conv3d_1/BiasAdd/ReadVariableOpReadVariableOp5functional_1_conv3d_1_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02.
,functional_1/conv3d_1/BiasAdd/ReadVariableOp?
functional_1/conv3d_1/BiasAddBiasAdd%functional_1/conv3d_1/Conv3D:output:04functional_1/conv3d_1/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2
functional_1/conv3d_1/BiasAdd?
functional_1/conv3d_1/EluElu&functional_1/conv3d_1/BiasAdd:output:0*
T0*3
_output_shapes!
:?????????2
functional_1/conv3d_1/Elu?
Afunctional_1/batch_normalization_1/moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2C
Afunctional_1/batch_normalization_1/moments/mean/reduction_indices?
/functional_1/batch_normalization_1/moments/meanMean'functional_1/conv3d_1/Elu:activations:0Jfunctional_1/batch_normalization_1/moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(21
/functional_1/batch_normalization_1/moments/mean?
7functional_1/batch_normalization_1/moments/StopGradientStopGradient8functional_1/batch_normalization_1/moments/mean:output:0*
T0**
_output_shapes
:29
7functional_1/batch_normalization_1/moments/StopGradient?
<functional_1/batch_normalization_1/moments/SquaredDifferenceSquaredDifference'functional_1/conv3d_1/Elu:activations:0@functional_1/batch_normalization_1/moments/StopGradient:output:0*
T0*3
_output_shapes!
:?????????2>
<functional_1/batch_normalization_1/moments/SquaredDifference?
Efunctional_1/batch_normalization_1/moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2G
Efunctional_1/batch_normalization_1/moments/variance/reduction_indices?
3functional_1/batch_normalization_1/moments/varianceMean@functional_1/batch_normalization_1/moments/SquaredDifference:z:0Nfunctional_1/batch_normalization_1/moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(25
3functional_1/batch_normalization_1/moments/variance?
2functional_1/batch_normalization_1/moments/SqueezeSqueeze8functional_1/batch_normalization_1/moments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 24
2functional_1/batch_normalization_1/moments/Squeeze?
4functional_1/batch_normalization_1/moments/Squeeze_1Squeeze<functional_1/batch_normalization_1/moments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 26
4functional_1/batch_normalization_1/moments/Squeeze_1?
8functional_1/batch_normalization_1/AssignMovingAvg/decayConst*M
_classC
A?loc:@functional_1/batch_normalization_1/AssignMovingAvg/1225960*
_output_shapes
: *
dtype0*
valueB
 *
?#<2:
8functional_1/batch_normalization_1/AssignMovingAvg/decay?
Afunctional_1/batch_normalization_1/AssignMovingAvg/ReadVariableOpReadVariableOp:functional_1_batch_normalization_1_assignmovingavg_1225960*
_output_shapes
:*
dtype02C
Afunctional_1/batch_normalization_1/AssignMovingAvg/ReadVariableOp?
6functional_1/batch_normalization_1/AssignMovingAvg/subSubIfunctional_1/batch_normalization_1/AssignMovingAvg/ReadVariableOp:value:0;functional_1/batch_normalization_1/moments/Squeeze:output:0*
T0*M
_classC
A?loc:@functional_1/batch_normalization_1/AssignMovingAvg/1225960*
_output_shapes
:28
6functional_1/batch_normalization_1/AssignMovingAvg/sub?
6functional_1/batch_normalization_1/AssignMovingAvg/mulMul:functional_1/batch_normalization_1/AssignMovingAvg/sub:z:0Afunctional_1/batch_normalization_1/AssignMovingAvg/decay:output:0*
T0*M
_classC
A?loc:@functional_1/batch_normalization_1/AssignMovingAvg/1225960*
_output_shapes
:28
6functional_1/batch_normalization_1/AssignMovingAvg/mul?
Ffunctional_1/batch_normalization_1/AssignMovingAvg/AssignSubVariableOpAssignSubVariableOp:functional_1_batch_normalization_1_assignmovingavg_1225960:functional_1/batch_normalization_1/AssignMovingAvg/mul:z:0B^functional_1/batch_normalization_1/AssignMovingAvg/ReadVariableOp*M
_classC
A?loc:@functional_1/batch_normalization_1/AssignMovingAvg/1225960*
_output_shapes
 *
dtype02H
Ffunctional_1/batch_normalization_1/AssignMovingAvg/AssignSubVariableOp?
:functional_1/batch_normalization_1/AssignMovingAvg_1/decayConst*O
_classE
CAloc:@functional_1/batch_normalization_1/AssignMovingAvg_1/1225966*
_output_shapes
: *
dtype0*
valueB
 *
?#<2<
:functional_1/batch_normalization_1/AssignMovingAvg_1/decay?
Cfunctional_1/batch_normalization_1/AssignMovingAvg_1/ReadVariableOpReadVariableOp<functional_1_batch_normalization_1_assignmovingavg_1_1225966*
_output_shapes
:*
dtype02E
Cfunctional_1/batch_normalization_1/AssignMovingAvg_1/ReadVariableOp?
8functional_1/batch_normalization_1/AssignMovingAvg_1/subSubKfunctional_1/batch_normalization_1/AssignMovingAvg_1/ReadVariableOp:value:0=functional_1/batch_normalization_1/moments/Squeeze_1:output:0*
T0*O
_classE
CAloc:@functional_1/batch_normalization_1/AssignMovingAvg_1/1225966*
_output_shapes
:2:
8functional_1/batch_normalization_1/AssignMovingAvg_1/sub?
8functional_1/batch_normalization_1/AssignMovingAvg_1/mulMul<functional_1/batch_normalization_1/AssignMovingAvg_1/sub:z:0Cfunctional_1/batch_normalization_1/AssignMovingAvg_1/decay:output:0*
T0*O
_classE
CAloc:@functional_1/batch_normalization_1/AssignMovingAvg_1/1225966*
_output_shapes
:2:
8functional_1/batch_normalization_1/AssignMovingAvg_1/mul?
Hfunctional_1/batch_normalization_1/AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOp<functional_1_batch_normalization_1_assignmovingavg_1_1225966<functional_1/batch_normalization_1/AssignMovingAvg_1/mul:z:0D^functional_1/batch_normalization_1/AssignMovingAvg_1/ReadVariableOp*O
_classE
CAloc:@functional_1/batch_normalization_1/AssignMovingAvg_1/1225966*
_output_shapes
 *
dtype02J
Hfunctional_1/batch_normalization_1/AssignMovingAvg_1/AssignSubVariableOp?
2functional_1/batch_normalization_1/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:24
2functional_1/batch_normalization_1/batchnorm/add/y?
0functional_1/batch_normalization_1/batchnorm/addAddV2=functional_1/batch_normalization_1/moments/Squeeze_1:output:0;functional_1/batch_normalization_1/batchnorm/add/y:output:0*
T0*
_output_shapes
:22
0functional_1/batch_normalization_1/batchnorm/add?
2functional_1/batch_normalization_1/batchnorm/RsqrtRsqrt4functional_1/batch_normalization_1/batchnorm/add:z:0*
T0*
_output_shapes
:24
2functional_1/batch_normalization_1/batchnorm/Rsqrt?
?functional_1/batch_normalization_1/batchnorm/mul/ReadVariableOpReadVariableOpHfunctional_1_batch_normalization_1_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02A
?functional_1/batch_normalization_1/batchnorm/mul/ReadVariableOp?
0functional_1/batch_normalization_1/batchnorm/mulMul6functional_1/batch_normalization_1/batchnorm/Rsqrt:y:0Gfunctional_1/batch_normalization_1/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:22
0functional_1/batch_normalization_1/batchnorm/mul?
2functional_1/batch_normalization_1/batchnorm/mul_1Mul'functional_1/conv3d_1/Elu:activations:04functional_1/batch_normalization_1/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????24
2functional_1/batch_normalization_1/batchnorm/mul_1?
2functional_1/batch_normalization_1/batchnorm/mul_2Mul;functional_1/batch_normalization_1/moments/Squeeze:output:04functional_1/batch_normalization_1/batchnorm/mul:z:0*
T0*
_output_shapes
:24
2functional_1/batch_normalization_1/batchnorm/mul_2?
;functional_1/batch_normalization_1/batchnorm/ReadVariableOpReadVariableOpDfunctional_1_batch_normalization_1_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02=
;functional_1/batch_normalization_1/batchnorm/ReadVariableOp?
0functional_1/batch_normalization_1/batchnorm/subSubCfunctional_1/batch_normalization_1/batchnorm/ReadVariableOp:value:06functional_1/batch_normalization_1/batchnorm/mul_2:z:0*
T0*
_output_shapes
:22
0functional_1/batch_normalization_1/batchnorm/sub?
2functional_1/batch_normalization_1/batchnorm/add_1AddV26functional_1/batch_normalization_1/batchnorm/mul_1:z:04functional_1/batch_normalization_1/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????24
2functional_1/batch_normalization_1/batchnorm/add_1?
+functional_1/conv3d_2/Conv3D/ReadVariableOpReadVariableOp4functional_1_conv3d_2_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02-
+functional_1/conv3d_2/Conv3D/ReadVariableOp?
functional_1/conv3d_2/Conv3DConv3D6functional_1/batch_normalization_1/batchnorm/add_1:z:03functional_1/conv3d_2/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
functional_1/conv3d_2/Conv3D?
,functional_1/conv3d_2/BiasAdd/ReadVariableOpReadVariableOp5functional_1_conv3d_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02.
,functional_1/conv3d_2/BiasAdd/ReadVariableOp?
functional_1/conv3d_2/BiasAddBiasAdd%functional_1/conv3d_2/Conv3D:output:04functional_1/conv3d_2/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2
functional_1/conv3d_2/BiasAdd?
functional_1/conv3d_2/EluElu&functional_1/conv3d_2/BiasAdd:output:0*
T0*3
_output_shapes!
:?????????2
functional_1/conv3d_2/Elu?
Afunctional_1/batch_normalization_2/moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2C
Afunctional_1/batch_normalization_2/moments/mean/reduction_indices?
/functional_1/batch_normalization_2/moments/meanMean'functional_1/conv3d_2/Elu:activations:0Jfunctional_1/batch_normalization_2/moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(21
/functional_1/batch_normalization_2/moments/mean?
7functional_1/batch_normalization_2/moments/StopGradientStopGradient8functional_1/batch_normalization_2/moments/mean:output:0*
T0**
_output_shapes
:29
7functional_1/batch_normalization_2/moments/StopGradient?
<functional_1/batch_normalization_2/moments/SquaredDifferenceSquaredDifference'functional_1/conv3d_2/Elu:activations:0@functional_1/batch_normalization_2/moments/StopGradient:output:0*
T0*3
_output_shapes!
:?????????2>
<functional_1/batch_normalization_2/moments/SquaredDifference?
Efunctional_1/batch_normalization_2/moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2G
Efunctional_1/batch_normalization_2/moments/variance/reduction_indices?
3functional_1/batch_normalization_2/moments/varianceMean@functional_1/batch_normalization_2/moments/SquaredDifference:z:0Nfunctional_1/batch_normalization_2/moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(25
3functional_1/batch_normalization_2/moments/variance?
2functional_1/batch_normalization_2/moments/SqueezeSqueeze8functional_1/batch_normalization_2/moments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 24
2functional_1/batch_normalization_2/moments/Squeeze?
4functional_1/batch_normalization_2/moments/Squeeze_1Squeeze<functional_1/batch_normalization_2/moments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 26
4functional_1/batch_normalization_2/moments/Squeeze_1?
8functional_1/batch_normalization_2/AssignMovingAvg/decayConst*M
_classC
A?loc:@functional_1/batch_normalization_2/AssignMovingAvg/1225999*
_output_shapes
: *
dtype0*
valueB
 *
?#<2:
8functional_1/batch_normalization_2/AssignMovingAvg/decay?
Afunctional_1/batch_normalization_2/AssignMovingAvg/ReadVariableOpReadVariableOp:functional_1_batch_normalization_2_assignmovingavg_1225999*
_output_shapes
:*
dtype02C
Afunctional_1/batch_normalization_2/AssignMovingAvg/ReadVariableOp?
6functional_1/batch_normalization_2/AssignMovingAvg/subSubIfunctional_1/batch_normalization_2/AssignMovingAvg/ReadVariableOp:value:0;functional_1/batch_normalization_2/moments/Squeeze:output:0*
T0*M
_classC
A?loc:@functional_1/batch_normalization_2/AssignMovingAvg/1225999*
_output_shapes
:28
6functional_1/batch_normalization_2/AssignMovingAvg/sub?
6functional_1/batch_normalization_2/AssignMovingAvg/mulMul:functional_1/batch_normalization_2/AssignMovingAvg/sub:z:0Afunctional_1/batch_normalization_2/AssignMovingAvg/decay:output:0*
T0*M
_classC
A?loc:@functional_1/batch_normalization_2/AssignMovingAvg/1225999*
_output_shapes
:28
6functional_1/batch_normalization_2/AssignMovingAvg/mul?
Ffunctional_1/batch_normalization_2/AssignMovingAvg/AssignSubVariableOpAssignSubVariableOp:functional_1_batch_normalization_2_assignmovingavg_1225999:functional_1/batch_normalization_2/AssignMovingAvg/mul:z:0B^functional_1/batch_normalization_2/AssignMovingAvg/ReadVariableOp*M
_classC
A?loc:@functional_1/batch_normalization_2/AssignMovingAvg/1225999*
_output_shapes
 *
dtype02H
Ffunctional_1/batch_normalization_2/AssignMovingAvg/AssignSubVariableOp?
:functional_1/batch_normalization_2/AssignMovingAvg_1/decayConst*O
_classE
CAloc:@functional_1/batch_normalization_2/AssignMovingAvg_1/1226005*
_output_shapes
: *
dtype0*
valueB
 *
?#<2<
:functional_1/batch_normalization_2/AssignMovingAvg_1/decay?
Cfunctional_1/batch_normalization_2/AssignMovingAvg_1/ReadVariableOpReadVariableOp<functional_1_batch_normalization_2_assignmovingavg_1_1226005*
_output_shapes
:*
dtype02E
Cfunctional_1/batch_normalization_2/AssignMovingAvg_1/ReadVariableOp?
8functional_1/batch_normalization_2/AssignMovingAvg_1/subSubKfunctional_1/batch_normalization_2/AssignMovingAvg_1/ReadVariableOp:value:0=functional_1/batch_normalization_2/moments/Squeeze_1:output:0*
T0*O
_classE
CAloc:@functional_1/batch_normalization_2/AssignMovingAvg_1/1226005*
_output_shapes
:2:
8functional_1/batch_normalization_2/AssignMovingAvg_1/sub?
8functional_1/batch_normalization_2/AssignMovingAvg_1/mulMul<functional_1/batch_normalization_2/AssignMovingAvg_1/sub:z:0Cfunctional_1/batch_normalization_2/AssignMovingAvg_1/decay:output:0*
T0*O
_classE
CAloc:@functional_1/batch_normalization_2/AssignMovingAvg_1/1226005*
_output_shapes
:2:
8functional_1/batch_normalization_2/AssignMovingAvg_1/mul?
Hfunctional_1/batch_normalization_2/AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOp<functional_1_batch_normalization_2_assignmovingavg_1_1226005<functional_1/batch_normalization_2/AssignMovingAvg_1/mul:z:0D^functional_1/batch_normalization_2/AssignMovingAvg_1/ReadVariableOp*O
_classE
CAloc:@functional_1/batch_normalization_2/AssignMovingAvg_1/1226005*
_output_shapes
 *
dtype02J
Hfunctional_1/batch_normalization_2/AssignMovingAvg_1/AssignSubVariableOp?
2functional_1/batch_normalization_2/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:24
2functional_1/batch_normalization_2/batchnorm/add/y?
0functional_1/batch_normalization_2/batchnorm/addAddV2=functional_1/batch_normalization_2/moments/Squeeze_1:output:0;functional_1/batch_normalization_2/batchnorm/add/y:output:0*
T0*
_output_shapes
:22
0functional_1/batch_normalization_2/batchnorm/add?
2functional_1/batch_normalization_2/batchnorm/RsqrtRsqrt4functional_1/batch_normalization_2/batchnorm/add:z:0*
T0*
_output_shapes
:24
2functional_1/batch_normalization_2/batchnorm/Rsqrt?
?functional_1/batch_normalization_2/batchnorm/mul/ReadVariableOpReadVariableOpHfunctional_1_batch_normalization_2_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02A
?functional_1/batch_normalization_2/batchnorm/mul/ReadVariableOp?
0functional_1/batch_normalization_2/batchnorm/mulMul6functional_1/batch_normalization_2/batchnorm/Rsqrt:y:0Gfunctional_1/batch_normalization_2/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:22
0functional_1/batch_normalization_2/batchnorm/mul?
2functional_1/batch_normalization_2/batchnorm/mul_1Mul'functional_1/conv3d_2/Elu:activations:04functional_1/batch_normalization_2/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????24
2functional_1/batch_normalization_2/batchnorm/mul_1?
2functional_1/batch_normalization_2/batchnorm/mul_2Mul;functional_1/batch_normalization_2/moments/Squeeze:output:04functional_1/batch_normalization_2/batchnorm/mul:z:0*
T0*
_output_shapes
:24
2functional_1/batch_normalization_2/batchnorm/mul_2?
;functional_1/batch_normalization_2/batchnorm/ReadVariableOpReadVariableOpDfunctional_1_batch_normalization_2_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02=
;functional_1/batch_normalization_2/batchnorm/ReadVariableOp?
0functional_1/batch_normalization_2/batchnorm/subSubCfunctional_1/batch_normalization_2/batchnorm/ReadVariableOp:value:06functional_1/batch_normalization_2/batchnorm/mul_2:z:0*
T0*
_output_shapes
:22
0functional_1/batch_normalization_2/batchnorm/sub?
2functional_1/batch_normalization_2/batchnorm/add_1AddV26functional_1/batch_normalization_2/batchnorm/mul_1:z:04functional_1/batch_normalization_2/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????24
2functional_1/batch_normalization_2/batchnorm/add_1?
+functional_1/conv3d_3/Conv3D/ReadVariableOpReadVariableOp4functional_1_conv3d_3_conv3d_readvariableop_resource**
_output_shapes
:*
dtype02-
+functional_1/conv3d_3/Conv3D/ReadVariableOp?
functional_1/conv3d_3/Conv3DConv3D6functional_1/batch_normalization_2/batchnorm/add_1:z:03functional_1/conv3d_3/Conv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
functional_1/conv3d_3/Conv3D?
,functional_1/conv3d_3/BiasAdd/ReadVariableOpReadVariableOp5functional_1_conv3d_3_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02.
,functional_1/conv3d_3/BiasAdd/ReadVariableOp?
functional_1/conv3d_3/BiasAddBiasAdd%functional_1/conv3d_3/Conv3D:output:04functional_1/conv3d_3/BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2
functional_1/conv3d_3/BiasAdd?
functional_1/conv3d_3/EluElu&functional_1/conv3d_3/BiasAdd:output:0*
T0*3
_output_shapes!
:?????????2
functional_1/conv3d_3/Elu?
Afunctional_1/batch_normalization_3/moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2C
Afunctional_1/batch_normalization_3/moments/mean/reduction_indices?
/functional_1/batch_normalization_3/moments/meanMean'functional_1/conv3d_3/Elu:activations:0Jfunctional_1/batch_normalization_3/moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(21
/functional_1/batch_normalization_3/moments/mean?
7functional_1/batch_normalization_3/moments/StopGradientStopGradient8functional_1/batch_normalization_3/moments/mean:output:0*
T0**
_output_shapes
:29
7functional_1/batch_normalization_3/moments/StopGradient?
<functional_1/batch_normalization_3/moments/SquaredDifferenceSquaredDifference'functional_1/conv3d_3/Elu:activations:0@functional_1/batch_normalization_3/moments/StopGradient:output:0*
T0*3
_output_shapes!
:?????????2>
<functional_1/batch_normalization_3/moments/SquaredDifference?
Efunctional_1/batch_normalization_3/moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2G
Efunctional_1/batch_normalization_3/moments/variance/reduction_indices?
3functional_1/batch_normalization_3/moments/varianceMean@functional_1/batch_normalization_3/moments/SquaredDifference:z:0Nfunctional_1/batch_normalization_3/moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(25
3functional_1/batch_normalization_3/moments/variance?
2functional_1/batch_normalization_3/moments/SqueezeSqueeze8functional_1/batch_normalization_3/moments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 24
2functional_1/batch_normalization_3/moments/Squeeze?
4functional_1/batch_normalization_3/moments/Squeeze_1Squeeze<functional_1/batch_normalization_3/moments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 26
4functional_1/batch_normalization_3/moments/Squeeze_1?
8functional_1/batch_normalization_3/AssignMovingAvg/decayConst*M
_classC
A?loc:@functional_1/batch_normalization_3/AssignMovingAvg/1226038*
_output_shapes
: *
dtype0*
valueB
 *
?#<2:
8functional_1/batch_normalization_3/AssignMovingAvg/decay?
Afunctional_1/batch_normalization_3/AssignMovingAvg/ReadVariableOpReadVariableOp:functional_1_batch_normalization_3_assignmovingavg_1226038*
_output_shapes
:*
dtype02C
Afunctional_1/batch_normalization_3/AssignMovingAvg/ReadVariableOp?
6functional_1/batch_normalization_3/AssignMovingAvg/subSubIfunctional_1/batch_normalization_3/AssignMovingAvg/ReadVariableOp:value:0;functional_1/batch_normalization_3/moments/Squeeze:output:0*
T0*M
_classC
A?loc:@functional_1/batch_normalization_3/AssignMovingAvg/1226038*
_output_shapes
:28
6functional_1/batch_normalization_3/AssignMovingAvg/sub?
6functional_1/batch_normalization_3/AssignMovingAvg/mulMul:functional_1/batch_normalization_3/AssignMovingAvg/sub:z:0Afunctional_1/batch_normalization_3/AssignMovingAvg/decay:output:0*
T0*M
_classC
A?loc:@functional_1/batch_normalization_3/AssignMovingAvg/1226038*
_output_shapes
:28
6functional_1/batch_normalization_3/AssignMovingAvg/mul?
Ffunctional_1/batch_normalization_3/AssignMovingAvg/AssignSubVariableOpAssignSubVariableOp:functional_1_batch_normalization_3_assignmovingavg_1226038:functional_1/batch_normalization_3/AssignMovingAvg/mul:z:0B^functional_1/batch_normalization_3/AssignMovingAvg/ReadVariableOp*M
_classC
A?loc:@functional_1/batch_normalization_3/AssignMovingAvg/1226038*
_output_shapes
 *
dtype02H
Ffunctional_1/batch_normalization_3/AssignMovingAvg/AssignSubVariableOp?
:functional_1/batch_normalization_3/AssignMovingAvg_1/decayConst*O
_classE
CAloc:@functional_1/batch_normalization_3/AssignMovingAvg_1/1226044*
_output_shapes
: *
dtype0*
valueB
 *
?#<2<
:functional_1/batch_normalization_3/AssignMovingAvg_1/decay?
Cfunctional_1/batch_normalization_3/AssignMovingAvg_1/ReadVariableOpReadVariableOp<functional_1_batch_normalization_3_assignmovingavg_1_1226044*
_output_shapes
:*
dtype02E
Cfunctional_1/batch_normalization_3/AssignMovingAvg_1/ReadVariableOp?
8functional_1/batch_normalization_3/AssignMovingAvg_1/subSubKfunctional_1/batch_normalization_3/AssignMovingAvg_1/ReadVariableOp:value:0=functional_1/batch_normalization_3/moments/Squeeze_1:output:0*
T0*O
_classE
CAloc:@functional_1/batch_normalization_3/AssignMovingAvg_1/1226044*
_output_shapes
:2:
8functional_1/batch_normalization_3/AssignMovingAvg_1/sub?
8functional_1/batch_normalization_3/AssignMovingAvg_1/mulMul<functional_1/batch_normalization_3/AssignMovingAvg_1/sub:z:0Cfunctional_1/batch_normalization_3/AssignMovingAvg_1/decay:output:0*
T0*O
_classE
CAloc:@functional_1/batch_normalization_3/AssignMovingAvg_1/1226044*
_output_shapes
:2:
8functional_1/batch_normalization_3/AssignMovingAvg_1/mul?
Hfunctional_1/batch_normalization_3/AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOp<functional_1_batch_normalization_3_assignmovingavg_1_1226044<functional_1/batch_normalization_3/AssignMovingAvg_1/mul:z:0D^functional_1/batch_normalization_3/AssignMovingAvg_1/ReadVariableOp*O
_classE
CAloc:@functional_1/batch_normalization_3/AssignMovingAvg_1/1226044*
_output_shapes
 *
dtype02J
Hfunctional_1/batch_normalization_3/AssignMovingAvg_1/AssignSubVariableOp?
2functional_1/batch_normalization_3/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:24
2functional_1/batch_normalization_3/batchnorm/add/y?
0functional_1/batch_normalization_3/batchnorm/addAddV2=functional_1/batch_normalization_3/moments/Squeeze_1:output:0;functional_1/batch_normalization_3/batchnorm/add/y:output:0*
T0*
_output_shapes
:22
0functional_1/batch_normalization_3/batchnorm/add?
2functional_1/batch_normalization_3/batchnorm/RsqrtRsqrt4functional_1/batch_normalization_3/batchnorm/add:z:0*
T0*
_output_shapes
:24
2functional_1/batch_normalization_3/batchnorm/Rsqrt?
?functional_1/batch_normalization_3/batchnorm/mul/ReadVariableOpReadVariableOpHfunctional_1_batch_normalization_3_batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02A
?functional_1/batch_normalization_3/batchnorm/mul/ReadVariableOp?
0functional_1/batch_normalization_3/batchnorm/mulMul6functional_1/batch_normalization_3/batchnorm/Rsqrt:y:0Gfunctional_1/batch_normalization_3/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:22
0functional_1/batch_normalization_3/batchnorm/mul?
2functional_1/batch_normalization_3/batchnorm/mul_1Mul'functional_1/conv3d_3/Elu:activations:04functional_1/batch_normalization_3/batchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????24
2functional_1/batch_normalization_3/batchnorm/mul_1?
2functional_1/batch_normalization_3/batchnorm/mul_2Mul;functional_1/batch_normalization_3/moments/Squeeze:output:04functional_1/batch_normalization_3/batchnorm/mul:z:0*
T0*
_output_shapes
:24
2functional_1/batch_normalization_3/batchnorm/mul_2?
;functional_1/batch_normalization_3/batchnorm/ReadVariableOpReadVariableOpDfunctional_1_batch_normalization_3_batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02=
;functional_1/batch_normalization_3/batchnorm/ReadVariableOp?
0functional_1/batch_normalization_3/batchnorm/subSubCfunctional_1/batch_normalization_3/batchnorm/ReadVariableOp:value:06functional_1/batch_normalization_3/batchnorm/mul_2:z:0*
T0*
_output_shapes
:22
0functional_1/batch_normalization_3/batchnorm/sub?
2functional_1/batch_normalization_3/batchnorm/add_1AddV26functional_1/batch_normalization_3/batchnorm/mul_1:z:04functional_1/batch_normalization_3/batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????24
2functional_1/batch_normalization_3/batchnorm/add_1?
(functional_1/average_pooling3d/AvgPool3D	AvgPool3D6functional_1/batch_normalization_3/batchnorm/add_1:z:0*
T0*3
_output_shapes!
:?????????*
ksize	
*
paddingVALID*
strides	
2*
(functional_1/average_pooling3d/AvgPool3D?
functional_1/flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"????   2
functional_1/flatten/Const?
functional_1/flatten/ReshapeReshape1functional_1/average_pooling3d/AvgPool3D:output:0#functional_1/flatten/Const:output:0*
T0*(
_output_shapes
:??????????
2
functional_1/flatten/Reshape?
"functional_1/dropout/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *UU??2$
"functional_1/dropout/dropout/Const?
 functional_1/dropout/dropout/MulMul%functional_1/flatten/Reshape:output:0+functional_1/dropout/dropout/Const:output:0*
T0*(
_output_shapes
:??????????
2"
 functional_1/dropout/dropout/Mul?
"functional_1/dropout/dropout/ShapeShape%functional_1/flatten/Reshape:output:0*
T0*
_output_shapes
:2$
"functional_1/dropout/dropout/Shape?
9functional_1/dropout/dropout/random_uniform/RandomUniformRandomUniform+functional_1/dropout/dropout/Shape:output:0*
T0*(
_output_shapes
:??????????
*
dtype02;
9functional_1/dropout/dropout/random_uniform/RandomUniform?
+functional_1/dropout/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2-
+functional_1/dropout/dropout/GreaterEqual/y?
)functional_1/dropout/dropout/GreaterEqualGreaterEqualBfunctional_1/dropout/dropout/random_uniform/RandomUniform:output:04functional_1/dropout/dropout/GreaterEqual/y:output:0*
T0*(
_output_shapes
:??????????
2+
)functional_1/dropout/dropout/GreaterEqual?
!functional_1/dropout/dropout/CastCast-functional_1/dropout/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*(
_output_shapes
:??????????
2#
!functional_1/dropout/dropout/Cast?
"functional_1/dropout/dropout/Mul_1Mul$functional_1/dropout/dropout/Mul:z:0%functional_1/dropout/dropout/Cast:y:0*
T0*(
_output_shapes
:??????????
2$
"functional_1/dropout/dropout/Mul_1?
)functional_1/layer1/MatMul/ReadVariableOpReadVariableOp2functional_1_layer1_matmul_readvariableop_resource* 
_output_shapes
:
?
?*
dtype02+
)functional_1/layer1/MatMul/ReadVariableOp?
functional_1/layer1/MatMulMatMul&functional_1/dropout/dropout/Mul_1:z:01functional_1/layer1/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2
functional_1/layer1/MatMul?
*functional_1/layer1/BiasAdd/ReadVariableOpReadVariableOp3functional_1_layer1_biasadd_readvariableop_resource*
_output_shapes	
:?*
dtype02,
*functional_1/layer1/BiasAdd/ReadVariableOp?
functional_1/layer1/BiasAddBiasAdd$functional_1/layer1/MatMul:product:02functional_1/layer1/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2
functional_1/layer1/BiasAdd?
functional_1/layer1/EluElu$functional_1/layer1/BiasAdd:output:0*
T0*(
_output_shapes
:??????????2
functional_1/layer1/Elu?
dense/MatMul/ReadVariableOpReadVariableOp$dense_matmul_readvariableop_resource*
_output_shapes
:	?d*
dtype02
dense/MatMul/ReadVariableOp?
dense/MatMulMatMul%functional_1/layer1/Elu:activations:0#dense/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????d2
dense/MatMul?
dense/BiasAdd/ReadVariableOpReadVariableOp%dense_biasadd_readvariableop_resource*
_output_shapes
:d*
dtype02
dense/BiasAdd/ReadVariableOp?
dense/BiasAddBiasAdddense/MatMul:product:0$dense/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????d2
dense/BiasAddg
	dense/EluEludense/BiasAdd:output:0*
T0*'
_output_shapes
:?????????d2
	dense/Elu?
dense_1/MatMul/ReadVariableOpReadVariableOp&dense_1_matmul_readvariableop_resource*
_output_shapes

:d
*
dtype02
dense_1/MatMul/ReadVariableOp?
dense_1/MatMulMatMuldense/Elu:activations:0%dense_1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????
2
dense_1/MatMul?
dense_1/BiasAdd/ReadVariableOpReadVariableOp'dense_1_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype02 
dense_1/BiasAdd/ReadVariableOp?
dense_1/BiasAddBiasAdddense_1/MatMul:product:0&dense_1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????
2
dense_1/BiasAddm
dense_1/EluEludense_1/BiasAdd:output:0*
T0*'
_output_shapes
:?????????
2
dense_1/Elu?
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOp$dense_matmul_readvariableop_resource*
_output_shapes
:	?d*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOp?
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?d2!
dense/kernel/Regularizer/Square?
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const?
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/Sum?
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2 
dense/kernel/Regularizer/mul/x?
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/mul?
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_1_matmul_readvariableop_resource*
_output_shapes

:d
*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp?
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d
2#
!dense_1/kernel/Regularizer/Square?
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const?
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/Sum?
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2"
 dense_1/kernel/Regularizer/mul/x?
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/mul?
IdentityIdentitydense_1/Elu:activations:0E^functional_1/batch_normalization/AssignMovingAvg/AssignSubVariableOpG^functional_1/batch_normalization/AssignMovingAvg_1/AssignSubVariableOpG^functional_1/batch_normalization_1/AssignMovingAvg/AssignSubVariableOpI^functional_1/batch_normalization_1/AssignMovingAvg_1/AssignSubVariableOpG^functional_1/batch_normalization_2/AssignMovingAvg/AssignSubVariableOpI^functional_1/batch_normalization_2/AssignMovingAvg_1/AssignSubVariableOpG^functional_1/batch_normalization_3/AssignMovingAvg/AssignSubVariableOpI^functional_1/batch_normalization_3/AssignMovingAvg_1/AssignSubVariableOp*
T0*'
_output_shapes
:?????????
2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????::::::::::::::::::::::::::::::2?
Dfunctional_1/batch_normalization/AssignMovingAvg/AssignSubVariableOpDfunctional_1/batch_normalization/AssignMovingAvg/AssignSubVariableOp2?
Ffunctional_1/batch_normalization/AssignMovingAvg_1/AssignSubVariableOpFfunctional_1/batch_normalization/AssignMovingAvg_1/AssignSubVariableOp2?
Ffunctional_1/batch_normalization_1/AssignMovingAvg/AssignSubVariableOpFfunctional_1/batch_normalization_1/AssignMovingAvg/AssignSubVariableOp2?
Hfunctional_1/batch_normalization_1/AssignMovingAvg_1/AssignSubVariableOpHfunctional_1/batch_normalization_1/AssignMovingAvg_1/AssignSubVariableOp2?
Ffunctional_1/batch_normalization_2/AssignMovingAvg/AssignSubVariableOpFfunctional_1/batch_normalization_2/AssignMovingAvg/AssignSubVariableOp2?
Hfunctional_1/batch_normalization_2/AssignMovingAvg_1/AssignSubVariableOpHfunctional_1/batch_normalization_2/AssignMovingAvg_1/AssignSubVariableOp2?
Ffunctional_1/batch_normalization_3/AssignMovingAvg/AssignSubVariableOpFfunctional_1/batch_normalization_3/AssignMovingAvg/AssignSubVariableOp2?
Hfunctional_1/batch_normalization_3/AssignMovingAvg_1/AssignSubVariableOpHfunctional_1/batch_normalization_3/AssignMovingAvg_1/AssignSubVariableOp:\ X
4
_output_shapes"
 :??????????
 
_user_specified_nameinputs
?*
?
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_1222947

inputs
assignmovingavg_1222922
assignmovingavg_1_1222928)
%batchnorm_mul_readvariableop_resource%
!batchnorm_readvariableop_resource
identity??#AssignMovingAvg/AssignSubVariableOp?%AssignMovingAvg_1/AssignSubVariableOp?
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2 
moments/mean/reduction_indices?
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/mean?
moments/StopGradientStopGradientmoments/mean:output:0*
T0**
_output_shapes
:2
moments/StopGradient?
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*3
_output_shapes!
:?????????2
moments/SquaredDifference?
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2$
"moments/variance/reduction_indices?
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/variance?
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze?
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze_1?
AssignMovingAvg/decayConst**
_class 
loc:@AssignMovingAvg/1222922*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg/decay?
AssignMovingAvg/ReadVariableOpReadVariableOpassignmovingavg_1222922*
_output_shapes
:*
dtype02 
AssignMovingAvg/ReadVariableOp?
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0**
_class 
loc:@AssignMovingAvg/1222922*
_output_shapes
:2
AssignMovingAvg/sub?
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0**
_class 
loc:@AssignMovingAvg/1222922*
_output_shapes
:2
AssignMovingAvg/mul?
#AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1222922AssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp**
_class 
loc:@AssignMovingAvg/1222922*
_output_shapes
 *
dtype02%
#AssignMovingAvg/AssignSubVariableOp?
AssignMovingAvg_1/decayConst*,
_class"
 loc:@AssignMovingAvg_1/1222928*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg_1/decay?
 AssignMovingAvg_1/ReadVariableOpReadVariableOpassignmovingavg_1_1222928*
_output_shapes
:*
dtype02"
 AssignMovingAvg_1/ReadVariableOp?
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1222928*
_output_shapes
:2
AssignMovingAvg_1/sub?
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1222928*
_output_shapes
:2
AssignMovingAvg_1/mul?
%AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1_1222928AssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*,
_class"
 loc:@AssignMovingAvg_1/1222928*
_output_shapes
 *
dtype02'
%AssignMovingAvg_1/AssignSubVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/mul_1{
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp?
batchnorm/subSub batchnorm/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0$^AssignMovingAvg/AssignSubVariableOp&^AssignMovingAvg_1/AssignSubVariableOp*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????::::2J
#AssignMovingAvg/AssignSubVariableOp#AssignMovingAvg/AssignSubVariableOp2N
%AssignMovingAvg_1/AssignSubVariableOp%AssignMovingAvg_1/AssignSubVariableOp:[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?
?
5__inference_batch_normalization_layer_call_fn_1226998

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *N
_output_shapes<
::8????????????????????????????????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *Y
fTRR
P__inference_batch_normalization_layer_call_and_return_conditional_losses_12222882
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????::::22
StatefulPartitionedCallStatefulPartitionedCall:v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
?9
?
I__inference_functional_3_layer_call_and_return_conditional_losses_1224678
input_1
input_2
input_3
sequential_1224566
sequential_1224568
sequential_1224570
sequential_1224572
sequential_1224574
sequential_1224576
sequential_1224578
sequential_1224580
sequential_1224582
sequential_1224584
sequential_1224586
sequential_1224588
sequential_1224590
sequential_1224592
sequential_1224594
sequential_1224596
sequential_1224598
sequential_1224600
sequential_1224602
sequential_1224604
sequential_1224606
sequential_1224608
sequential_1224610
sequential_1224612
sequential_1224614
sequential_1224616
sequential_1224618
sequential_1224620
sequential_1224622
sequential_1224624
dense_2_1224660
dense_2_1224662
identity??dense_2/StatefulPartitionedCall?"sequential/StatefulPartitionedCall?$sequential/StatefulPartitionedCall_1?
"sequential/StatefulPartitionedCallStatefulPartitionedCallinput_1sequential_1224566sequential_1224568sequential_1224570sequential_1224572sequential_1224574sequential_1224576sequential_1224578sequential_1224580sequential_1224582sequential_1224584sequential_1224586sequential_1224588sequential_1224590sequential_1224592sequential_1224594sequential_1224596sequential_1224598sequential_1224600sequential_1224602sequential_1224604sequential_1224606sequential_1224608sequential_1224610sequential_1224612sequential_1224614sequential_1224616sequential_1224618sequential_1224620sequential_1224622sequential_1224624**
Tin#
!2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
*@
_read_only_resource_inputs"
 	
*0
config_proto 

CPU

GPU2*0J 8? *P
fKRI
G__inference_sequential_layer_call_and_return_conditional_losses_12241862$
"sequential/StatefulPartitionedCall?
$sequential/StatefulPartitionedCall_1StatefulPartitionedCallinput_2sequential_1224566sequential_1224568sequential_1224570sequential_1224572sequential_1224574sequential_1224576sequential_1224578sequential_1224580sequential_1224582sequential_1224584sequential_1224586sequential_1224588sequential_1224590sequential_1224592sequential_1224594sequential_1224596sequential_1224598sequential_1224600sequential_1224602sequential_1224604sequential_1224606sequential_1224608sequential_1224610sequential_1224612sequential_1224614sequential_1224616sequential_1224618sequential_1224620sequential_1224622sequential_1224624**
Tin#
!2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
*@
_read_only_resource_inputs"
 	
*0
config_proto 

CPU

GPU2*0J 8? *P
fKRI
G__inference_sequential_layer_call_and_return_conditional_losses_12241862&
$sequential/StatefulPartitionedCall_1?
lambda/PartitionedCallPartitionedCall+sequential/StatefulPartitionedCall:output:0-sequential/StatefulPartitionedCall_1:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *L
fGRE
C__inference_lambda_layer_call_and_return_conditional_losses_12244912
lambda/PartitionedCall?
concatenate/PartitionedCallPartitionedCalllambda/PartitionedCall:output:0input_3*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *Q
fLRJ
H__inference_concatenate_layer_call_and_return_conditional_losses_12245132
concatenate/PartitionedCall?
dense_2/StatefulPartitionedCallStatefulPartitionedCall$concatenate/PartitionedCall:output:0dense_2_1224660dense_2_1224662*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_2_layer_call_and_return_conditional_losses_12245322!
dense_2/StatefulPartitionedCall?
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_1224618*
_output_shapes
:	?d*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOp?
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?d2!
dense/kernel/Regularizer/Square?
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const?
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/Sum?
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2 
dense/kernel/Regularizer/mul/x?
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/mul?
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_1224622*
_output_shapes

:d
*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp?
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d
2#
!dense_1/kernel/Regularizer/Square?
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const?
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/Sum?
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2"
 dense_1/kernel/Regularizer/mul/x?
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/mul?
IdentityIdentity(dense_2/StatefulPartitionedCall:output:0 ^dense_2/StatefulPartitionedCall#^sequential/StatefulPartitionedCall%^sequential/StatefulPartitionedCall_1*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????:??????????:?????????
::::::::::::::::::::::::::::::::2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2H
"sequential/StatefulPartitionedCall"sequential/StatefulPartitionedCall2L
$sequential/StatefulPartitionedCall_1$sequential/StatefulPartitionedCall_1:] Y
4
_output_shapes"
 :??????????
!
_user_specified_name	input_1:]Y
4
_output_shapes"
 :??????????
!
_user_specified_name	input_2:PL
'
_output_shapes
:?????????

!
_user_specified_name	input_3
?
?
D__inference_dense_2_layer_call_and_return_conditional_losses_1226418

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity??
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
MatMul?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2	
BiasAddd
IdentityIdentityBiasAdd:output:0*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*.
_input_shapes
:?????????:::O K
'
_output_shapes
:?????????
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_1227169

inputs%
!batchnorm_readvariableop_resource)
%batchnorm_mul_readvariableop_resource'
#batchnorm_readvariableop_1_resource'
#batchnorm_readvariableop_2_resource
identity??
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/mul_1?
batchnorm/ReadVariableOp_1ReadVariableOp#batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_1?
batchnorm/mul_2Mul"batchnorm/ReadVariableOp_1:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOp_2ReadVariableOp#batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_2?
batchnorm/subSub"batchnorm/ReadVariableOp_2:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0*
T0*N
_output_shapes<
::8????????????????????????????????????2

Identity"
identityIdentity:output:0*]
_input_shapesL
J:8????????????????????????????????????:::::v r
N
_output_shapes<
::8????????????????????????????????????
 
_user_specified_nameinputs
?
|
'__inference_dense_layer_call_fn_1226856

inputs
unknown
	unknown_0
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????d*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *K
fFRD
B__inference_dense_layer_call_and_return_conditional_losses_12238192
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????d2

Identity"
identityIdentity:output:0*/
_input_shapes
:??????????::22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs
?
~
)__inference_dense_1_layer_call_fn_1226888

inputs
unknown
	unknown_0
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_1_layer_call_and_return_conditional_losses_12238522
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????
2

Identity"
identityIdentity:output:0*.
_input_shapes
:?????????d::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:?????????d
 
_user_specified_nameinputs
?
b
)__inference_dropout_layer_call_fn_1227678

inputs
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????
* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dropout_layer_call_and_return_conditional_losses_12232662
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*(
_output_shapes
:??????????
2

Identity"
identityIdentity:output:0*'
_input_shapes
:??????????
22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:??????????

 
_user_specified_nameinputs
?
?
7__inference_batch_normalization_3_layer_call_fn_1227632

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_12231832
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????::::22
StatefulPartitionedCallStatefulPartitionedCall:[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?
?
.__inference_functional_1_layer_call_fn_1226824

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16

unknown_17

unknown_18

unknown_19

unknown_20

unknown_21

unknown_22

unknown_23

unknown_24
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19
unknown_20
unknown_21
unknown_22
unknown_23
unknown_24*&
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*<
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8? *R
fMRK
I__inference_functional_1_layer_call_and_return_conditional_losses_12235762
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????::::::::::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:\ X
4
_output_shapes"
 :??????????
 
_user_specified_nameinputs
?*
?
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_1223183

inputs
assignmovingavg_1223158
assignmovingavg_1_1223164)
%batchnorm_mul_readvariableop_resource%
!batchnorm_readvariableop_resource
identity??#AssignMovingAvg/AssignSubVariableOp?%AssignMovingAvg_1/AssignSubVariableOp?
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2 
moments/mean/reduction_indices?
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/mean?
moments/StopGradientStopGradientmoments/mean:output:0*
T0**
_output_shapes
:2
moments/StopGradient?
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*3
_output_shapes!
:?????????2
moments/SquaredDifference?
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*%
valueB"             2$
"moments/variance/reduction_indices?
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0**
_output_shapes
:*
	keep_dims(2
moments/variance?
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze?
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 2
moments/Squeeze_1?
AssignMovingAvg/decayConst**
_class 
loc:@AssignMovingAvg/1223158*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg/decay?
AssignMovingAvg/ReadVariableOpReadVariableOpassignmovingavg_1223158*
_output_shapes
:*
dtype02 
AssignMovingAvg/ReadVariableOp?
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0**
_class 
loc:@AssignMovingAvg/1223158*
_output_shapes
:2
AssignMovingAvg/sub?
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0**
_class 
loc:@AssignMovingAvg/1223158*
_output_shapes
:2
AssignMovingAvg/mul?
#AssignMovingAvg/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1223158AssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp**
_class 
loc:@AssignMovingAvg/1223158*
_output_shapes
 *
dtype02%
#AssignMovingAvg/AssignSubVariableOp?
AssignMovingAvg_1/decayConst*,
_class"
 loc:@AssignMovingAvg_1/1223164*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg_1/decay?
 AssignMovingAvg_1/ReadVariableOpReadVariableOpassignmovingavg_1_1223164*
_output_shapes
:*
dtype02"
 AssignMovingAvg_1/ReadVariableOp?
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1223164*
_output_shapes
:2
AssignMovingAvg_1/sub?
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*,
_class"
 loc:@AssignMovingAvg_1/1223164*
_output_shapes
:2
AssignMovingAvg_1/mul?
%AssignMovingAvg_1/AssignSubVariableOpAssignSubVariableOpassignmovingavg_1_1223164AssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*,
_class"
 loc:@AssignMovingAvg_1/1223164*
_output_shapes
 *
dtype02'
%AssignMovingAvg_1/AssignSubVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/mul_1{
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp?
batchnorm/subSub batchnorm/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/add_1?
IdentityIdentitybatchnorm/add_1:z:0$^AssignMovingAvg/AssignSubVariableOp&^AssignMovingAvg_1/AssignSubVariableOp*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????::::2J
#AssignMovingAvg/AssignSubVariableOp#AssignMovingAvg/AssignSubVariableOp2N
%AssignMovingAvg_1/AssignSubVariableOp%AssignMovingAvg_1/AssignSubVariableOp:[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?
r
H__inference_concatenate_layer_call_and_return_conditional_losses_1224513

inputs
inputs_1
identity\
concat/axisConst*
_output_shapes
: *
dtype0*
value	B :2
concat/axis
concatConcatV2inputsinputs_1concat/axis:output:0*
N*
T0*'
_output_shapes
:?????????2
concatc
IdentityIdentityconcat:output:0*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*9
_input_shapes(
&:?????????
:?????????
:O K
'
_output_shapes
:?????????

 
_user_specified_nameinputs:OK
'
_output_shapes
:?????????

 
_user_specified_nameinputs
?B
?	
I__inference_functional_1_layer_call_and_return_conditional_losses_1223576

inputs
conv3d_1223511
conv3d_1223513
batch_normalization_1223516
batch_normalization_1223518
batch_normalization_1223520
batch_normalization_1223522
conv3d_1_1223525
conv3d_1_1223527!
batch_normalization_1_1223530!
batch_normalization_1_1223532!
batch_normalization_1_1223534!
batch_normalization_1_1223536
conv3d_2_1223539
conv3d_2_1223541!
batch_normalization_2_1223544!
batch_normalization_2_1223546!
batch_normalization_2_1223548!
batch_normalization_2_1223550
conv3d_3_1223553
conv3d_3_1223555!
batch_normalization_3_1223558!
batch_normalization_3_1223560!
batch_normalization_3_1223562!
batch_normalization_3_1223564
layer1_1223570
layer1_1223572
identity??+batch_normalization/StatefulPartitionedCall?-batch_normalization_1/StatefulPartitionedCall?-batch_normalization_2/StatefulPartitionedCall?-batch_normalization_3/StatefulPartitionedCall?conv3d/StatefulPartitionedCall? conv3d_1/StatefulPartitionedCall? conv3d_2/StatefulPartitionedCall? conv3d_3/StatefulPartitionedCall?layer1/StatefulPartitionedCall?
conv3d/StatefulPartitionedCallStatefulPartitionedCallinputsconv3d_1223511conv3d_1223513*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *L
fGRE
C__inference_conv3d_layer_call_and_return_conditional_losses_12227782 
conv3d/StatefulPartitionedCall?
+batch_normalization/StatefulPartitionedCallStatefulPartitionedCall'conv3d/StatefulPartitionedCall:output:0batch_normalization_1223516batch_normalization_1223518batch_normalization_1223520batch_normalization_1223522*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *Y
fTRR
P__inference_batch_normalization_layer_call_and_return_conditional_losses_12228492-
+batch_normalization/StatefulPartitionedCall?
 conv3d_1/StatefulPartitionedCallStatefulPartitionedCall4batch_normalization/StatefulPartitionedCall:output:0conv3d_1_1223525conv3d_1_1223527*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv3d_1_layer_call_and_return_conditional_losses_12228962"
 conv3d_1/StatefulPartitionedCall?
-batch_normalization_1/StatefulPartitionedCallStatefulPartitionedCall)conv3d_1/StatefulPartitionedCall:output:0batch_normalization_1_1223530batch_normalization_1_1223532batch_normalization_1_1223534batch_normalization_1_1223536*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_12229672/
-batch_normalization_1/StatefulPartitionedCall?
 conv3d_2/StatefulPartitionedCallStatefulPartitionedCall6batch_normalization_1/StatefulPartitionedCall:output:0conv3d_2_1223539conv3d_2_1223541*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv3d_2_layer_call_and_return_conditional_losses_12230142"
 conv3d_2/StatefulPartitionedCall?
-batch_normalization_2/StatefulPartitionedCallStatefulPartitionedCall)conv3d_2/StatefulPartitionedCall:output:0batch_normalization_2_1223544batch_normalization_2_1223546batch_normalization_2_1223548batch_normalization_2_1223550*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_12230852/
-batch_normalization_2/StatefulPartitionedCall?
 conv3d_3/StatefulPartitionedCallStatefulPartitionedCall6batch_normalization_2/StatefulPartitionedCall:output:0conv3d_3_1223553conv3d_3_1223555*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv3d_3_layer_call_and_return_conditional_losses_12231322"
 conv3d_3/StatefulPartitionedCall?
-batch_normalization_3/StatefulPartitionedCallStatefulPartitionedCall)conv3d_3/StatefulPartitionedCall:output:0batch_normalization_3_1223558batch_normalization_3_1223560batch_normalization_3_1223562batch_normalization_3_1223564*
Tin	
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_12232032/
-batch_normalization_3/StatefulPartitionedCall?
!average_pooling3d/PartitionedCallPartitionedCall6batch_normalization_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *3
_output_shapes!
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *W
fRRP
N__inference_average_pooling3d_layer_call_and_return_conditional_losses_12227582#
!average_pooling3d/PartitionedCall?
flatten/PartitionedCallPartitionedCall*average_pooling3d/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????
* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_flatten_layer_call_and_return_conditional_losses_12232462
flatten/PartitionedCall?
dropout/PartitionedCallPartitionedCall flatten/PartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????
* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dropout_layer_call_and_return_conditional_losses_12232712
dropout/PartitionedCall?
layer1/StatefulPartitionedCallStatefulPartitionedCall dropout/PartitionedCall:output:0layer1_1223570layer1_1223572*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *L
fGRE
C__inference_layer1_layer_call_and_return_conditional_losses_12232952 
layer1/StatefulPartitionedCall?
IdentityIdentity'layer1/StatefulPartitionedCall:output:0,^batch_normalization/StatefulPartitionedCall.^batch_normalization_1/StatefulPartitionedCall.^batch_normalization_2/StatefulPartitionedCall.^batch_normalization_3/StatefulPartitionedCall^conv3d/StatefulPartitionedCall!^conv3d_1/StatefulPartitionedCall!^conv3d_2/StatefulPartitionedCall!^conv3d_3/StatefulPartitionedCall^layer1/StatefulPartitionedCall*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????::::::::::::::::::::::::::2Z
+batch_normalization/StatefulPartitionedCall+batch_normalization/StatefulPartitionedCall2^
-batch_normalization_1/StatefulPartitionedCall-batch_normalization_1/StatefulPartitionedCall2^
-batch_normalization_2/StatefulPartitionedCall-batch_normalization_2/StatefulPartitionedCall2^
-batch_normalization_3/StatefulPartitionedCall-batch_normalization_3/StatefulPartitionedCall2@
conv3d/StatefulPartitionedCallconv3d/StatefulPartitionedCall2D
 conv3d_1/StatefulPartitionedCall conv3d_1/StatefulPartitionedCall2D
 conv3d_2/StatefulPartitionedCall conv3d_2/StatefulPartitionedCall2D
 conv3d_3/StatefulPartitionedCall conv3d_3/StatefulPartitionedCall2@
layer1/StatefulPartitionedCalllayer1/StatefulPartitionedCall:\ X
4
_output_shapes"
 :??????????
 
_user_specified_nameinputs
?
m
C__inference_lambda_layer_call_and_return_conditional_losses_1224491

inputs
inputs_1
identityU
subSubinputsinputs_1*
T0*'
_output_shapes
:?????????
2
subL
AbsAbssub:z:0*
T0*'
_output_shapes
:?????????
2
Abs[
IdentityIdentityAbs:y:0*
T0*'
_output_shapes
:?????????
2

Identity"
identityIdentity:output:0*9
_input_shapes(
&:?????????
:?????????
:O K
'
_output_shapes
:?????????

 
_user_specified_nameinputs:OK
'
_output_shapes
:?????????

 
_user_specified_nameinputs
?
?
,__inference_sequential_layer_call_fn_1224249
functional_1_input
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16

unknown_17

unknown_18

unknown_19

unknown_20

unknown_21

unknown_22

unknown_23

unknown_24

unknown_25

unknown_26

unknown_27

unknown_28
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallfunctional_1_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19
unknown_20
unknown_21
unknown_22
unknown_23
unknown_24
unknown_25
unknown_26
unknown_27
unknown_28**
Tin#
!2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
*@
_read_only_resource_inputs"
 	
*0
config_proto 

CPU

GPU2*0J 8? *P
fKRI
G__inference_sequential_layer_call_and_return_conditional_losses_12241862
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????
2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????::::::::::::::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:h d
4
_output_shapes"
 :??????????
,
_user_specified_namefunctional_1_input
?
?
.__inference_functional_3_layer_call_fn_1224867
input_1
input_2
input_3
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
	unknown_8
	unknown_9

unknown_10

unknown_11

unknown_12

unknown_13

unknown_14

unknown_15

unknown_16

unknown_17

unknown_18

unknown_19

unknown_20

unknown_21

unknown_22

unknown_23

unknown_24

unknown_25

unknown_26

unknown_27

unknown_28

unknown_29

unknown_30
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinput_1input_2input_3unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19
unknown_20
unknown_21
unknown_22
unknown_23
unknown_24
unknown_25
unknown_26
unknown_27
unknown_28
unknown_29
unknown_30*.
Tin'
%2#*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????*:
_read_only_resource_inputs
	
 !"*0
config_proto 

CPU

GPU2*0J 8? *R
fMRK
I__inference_functional_3_layer_call_and_return_conditional_losses_12248002
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????:??????????:?????????
::::::::::::::::::::::::::::::::22
StatefulPartitionedCallStatefulPartitionedCall:] Y
4
_output_shapes"
 :??????????
!
_user_specified_name	input_1:]Y
4
_output_shapes"
 :??????????
!
_user_specified_name	input_2:PL
'
_output_shapes
:?????????

!
_user_specified_name	input_3
?	
?
E__inference_conv3d_3_layer_call_and_return_conditional_losses_1223132

inputs"
conv3d_readvariableop_resource#
biasadd_readvariableop_resource
identity??
Conv3D/ReadVariableOpReadVariableOpconv3d_readvariableop_resource**
_output_shapes
:*
dtype02
Conv3D/ReadVariableOp?
Conv3DConv3DinputsConv3D/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????*
paddingVALID*
strides	
2
Conv3D?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddConv3D:output:0BiasAdd/ReadVariableOp:value:0*
T0*3
_output_shapes!
:?????????2	
BiasAdda
EluEluBiasAdd:output:0*
T0*3
_output_shapes!
:?????????2
Eluq
IdentityIdentityElu:activations:0*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*:
_input_shapes)
':?????????:::[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs
?9
?
I__inference_functional_3_layer_call_and_return_conditional_losses_1224800

inputs
inputs_1
inputs_2
sequential_1224688
sequential_1224690
sequential_1224692
sequential_1224694
sequential_1224696
sequential_1224698
sequential_1224700
sequential_1224702
sequential_1224704
sequential_1224706
sequential_1224708
sequential_1224710
sequential_1224712
sequential_1224714
sequential_1224716
sequential_1224718
sequential_1224720
sequential_1224722
sequential_1224724
sequential_1224726
sequential_1224728
sequential_1224730
sequential_1224732
sequential_1224734
sequential_1224736
sequential_1224738
sequential_1224740
sequential_1224742
sequential_1224744
sequential_1224746
dense_2_1224782
dense_2_1224784
identity??dense_2/StatefulPartitionedCall?"sequential/StatefulPartitionedCall?$sequential/StatefulPartitionedCall_1?
"sequential/StatefulPartitionedCallStatefulPartitionedCallinputssequential_1224688sequential_1224690sequential_1224692sequential_1224694sequential_1224696sequential_1224698sequential_1224700sequential_1224702sequential_1224704sequential_1224706sequential_1224708sequential_1224710sequential_1224712sequential_1224714sequential_1224716sequential_1224718sequential_1224720sequential_1224722sequential_1224724sequential_1224726sequential_1224728sequential_1224730sequential_1224732sequential_1224734sequential_1224736sequential_1224738sequential_1224740sequential_1224742sequential_1224744sequential_1224746**
Tin#
!2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
*8
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *P
fKRI
G__inference_sequential_layer_call_and_return_conditional_losses_12240422$
"sequential/StatefulPartitionedCall?
$sequential/StatefulPartitionedCall_1StatefulPartitionedCallinputs_1sequential_1224688sequential_1224690sequential_1224692sequential_1224694sequential_1224696sequential_1224698sequential_1224700sequential_1224702sequential_1224704sequential_1224706sequential_1224708sequential_1224710sequential_1224712sequential_1224714sequential_1224716sequential_1224718sequential_1224720sequential_1224722sequential_1224724sequential_1224726sequential_1224728sequential_1224730sequential_1224732sequential_1224734sequential_1224736sequential_1224738sequential_1224740sequential_1224742sequential_1224744sequential_1224746#^sequential/StatefulPartitionedCall**
Tin#
!2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
*8
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *P
fKRI
G__inference_sequential_layer_call_and_return_conditional_losses_12240422&
$sequential/StatefulPartitionedCall_1?
lambda/PartitionedCallPartitionedCall+sequential/StatefulPartitionedCall:output:0-sequential/StatefulPartitionedCall_1:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????
* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *L
fGRE
C__inference_lambda_layer_call_and_return_conditional_losses_12244842
lambda/PartitionedCall?
concatenate/PartitionedCallPartitionedCalllambda/PartitionedCall:output:0inputs_2*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *Q
fLRJ
H__inference_concatenate_layer_call_and_return_conditional_losses_12245132
concatenate/PartitionedCall?
dense_2/StatefulPartitionedCallStatefulPartitionedCall$concatenate/PartitionedCall:output:0dense_2_1224782dense_2_1224784*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_2_layer_call_and_return_conditional_losses_12245322!
dense_2/StatefulPartitionedCall?
.dense/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_1224740*
_output_shapes
:	?d*
dtype020
.dense/kernel/Regularizer/Square/ReadVariableOp?
dense/kernel/Regularizer/SquareSquare6dense/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?d2!
dense/kernel/Regularizer/Square?
dense/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2 
dense/kernel/Regularizer/Const?
dense/kernel/Regularizer/SumSum#dense/kernel/Regularizer/Square:y:0'dense/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/Sum?
dense/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2 
dense/kernel/Regularizer/mul/x?
dense/kernel/Regularizer/mulMul'dense/kernel/Regularizer/mul/x:output:0%dense/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense/kernel/Regularizer/mul?
0dense_1/kernel/Regularizer/Square/ReadVariableOpReadVariableOpsequential_1224744*
_output_shapes

:d
*
dtype022
0dense_1/kernel/Regularizer/Square/ReadVariableOp?
!dense_1/kernel/Regularizer/SquareSquare8dense_1/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes

:d
2#
!dense_1/kernel/Regularizer/Square?
 dense_1/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_1/kernel/Regularizer/Const?
dense_1/kernel/Regularizer/SumSum%dense_1/kernel/Regularizer/Square:y:0)dense_1/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/Sum?
 dense_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2"
 dense_1/kernel/Regularizer/mul/x?
dense_1/kernel/Regularizer/mulMul)dense_1/kernel/Regularizer/mul/x:output:0'dense_1/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_1/kernel/Regularizer/mul?
IdentityIdentity(dense_2/StatefulPartitionedCall:output:0 ^dense_2/StatefulPartitionedCall#^sequential/StatefulPartitionedCall%^sequential/StatefulPartitionedCall_1*
T0*'
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*?
_input_shapes?
?:??????????:??????????:?????????
::::::::::::::::::::::::::::::::2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2H
"sequential/StatefulPartitionedCall"sequential/StatefulPartitionedCall2L
$sequential/StatefulPartitionedCall_1$sequential/StatefulPartitionedCall_1:\ X
4
_output_shapes"
 :??????????
 
_user_specified_nameinputs:\X
4
_output_shapes"
 :??????????
 
_user_specified_nameinputs:OK
'
_output_shapes
:?????????

 
_user_specified_nameinputs
?
o
C__inference_lambda_layer_call_and_return_conditional_losses_1226383
inputs_0
inputs_1
identityW
subSubinputs_0inputs_1*
T0*'
_output_shapes
:?????????
2
subL
AbsAbssub:z:0*
T0*'
_output_shapes
:?????????
2
Abs[
IdentityIdentityAbs:y:0*
T0*'
_output_shapes
:?????????
2

Identity"
identityIdentity:output:0*9
_input_shapes(
&:?????????
:?????????
:Q M
'
_output_shapes
:?????????

"
_user_specified_name
inputs/0:QM
'
_output_shapes
:?????????

"
_user_specified_name
inputs/1
?
?
P__inference_batch_normalization_layer_call_and_return_conditional_losses_1227067

inputs%
!batchnorm_readvariableop_resource)
%batchnorm_mul_readvariableop_resource'
#batchnorm_readvariableop_1_resource'
#batchnorm_readvariableop_2_resource
identity??
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes
:2
batchnorm/addc
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes
:2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes
:*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes
:2
batchnorm/mul?
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/mul_1?
batchnorm/ReadVariableOp_1ReadVariableOp#batchnorm_readvariableop_1_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_1?
batchnorm/mul_2Mul"batchnorm/ReadVariableOp_1:value:0batchnorm/mul:z:0*
T0*
_output_shapes
:2
batchnorm/mul_2?
batchnorm/ReadVariableOp_2ReadVariableOp#batchnorm_readvariableop_2_resource*
_output_shapes
:*
dtype02
batchnorm/ReadVariableOp_2?
batchnorm/subSub"batchnorm/ReadVariableOp_2:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes
:2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*3
_output_shapes!
:?????????2
batchnorm/add_1s
IdentityIdentitybatchnorm/add_1:z:0*
T0*3
_output_shapes!
:?????????2

Identity"
identityIdentity:output:0*B
_input_shapes1
/:?????????:::::[ W
3
_output_shapes!
:?????????
 
_user_specified_nameinputs"?L
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*?
serving_default?
H
input_1=
serving_default_input_1:0??????????
H
input_2=
serving_default_input_2:0??????????
;
input_30
serving_default_input_3:0?????????
;
dense_20
StatefulPartitionedCall:0?????????tensorflow/serving/predict:??
??
layer-0
layer-1
layer_with_weights-0
layer-2
layer-3
layer-4
layer-5
layer_with_weights-1
layer-6
	optimizer
	trainable_variables

regularization_losses
	variables
	keras_api

signatures
+?&call_and_return_all_conditional_losses
?__call__
?_default_save_signature"ص
_tf_keras_network??{"class_name": "Functional", "name": "functional_3", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "must_restore_from_config": false, "config": {"name": "functional_3", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 24, 24, 24, 167]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_1"}, "name": "input_1", "inbound_nodes": []}, {"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 24, 24, 24, 167]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_2"}, "name": "input_2", "inbound_nodes": []}, {"class_name": "Sequential", "config": {"name": "sequential", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 24, 24, 24, 167]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "functional_1_input"}}, {"class_name": "Functional", "config": {"name": "functional_1", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 24, 24, 24, 167]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_1"}, "name": "input_1", "inbound_nodes": []}, {"class_name": "Conv3D", "config": {"name": "conv3d", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 24, 24, 24, 167]}, "dtype": "float32", "filters": 20, "kernel_size": {"class_name": "__tuple__", "items": [1, 1, 1]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d", "inbound_nodes": [[["input_1", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization", "inbound_nodes": [[["conv3d", 0, 0, {}]]]}, {"class_name": "Conv3D", "config": {"name": "conv3d_1", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 24, 24, 24, 20]}, "dtype": "float32", "filters": 20, "kernel_size": {"class_name": "__tuple__", "items": [3, 3, 3]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d_1", "inbound_nodes": [[["batch_normalization", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization_1", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization_1", "inbound_nodes": [[["conv3d_1", 0, 0, {}]]]}, {"class_name": "Conv3D", "config": {"name": "conv3d_2", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 22, 22, 22, 20]}, "dtype": "float32", "filters": 30, "kernel_size": {"class_name": "__tuple__", "items": [4, 4, 4]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d_2", "inbound_nodes": [[["batch_normalization_1", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization_2", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization_2", "inbound_nodes": [[["conv3d_2", 0, 0, {}]]]}, {"class_name": "Conv3D", "config": {"name": "conv3d_3", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 19, 19, 19, 30]}, "dtype": "float32", "filters": 20, "kernel_size": {"class_name": "__tuple__", "items": [4, 4, 4]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d_3", "inbound_nodes": [[["batch_normalization_2", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization_3", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization_3", "inbound_nodes": [[["conv3d_3", 0, 0, {}]]]}, {"class_name": "AveragePooling3D", "config": {"name": "average_pooling3d", "trainable": true, "dtype": "float32", "pool_size": {"class_name": "__tuple__", "items": [4, 4, 4]}, "padding": "valid", "strides": {"class_name": "__tuple__", "items": [4, 4, 4]}, "data_format": "channels_last"}, "name": "average_pooling3d", "inbound_nodes": [[["batch_normalization_3", 0, 0, {}]]]}, {"class_name": "Flatten", "config": {"name": "flatten", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "name": "flatten", "inbound_nodes": [[["average_pooling3d", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "dropout", "trainable": true, "dtype": "float32", "rate": 0.4, "noise_shape": null, "seed": null}, "name": "dropout", "inbound_nodes": [[["flatten", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "layer1", "trainable": true, "dtype": "float32", "units": 200, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": {"class_name": "MaxNorm", "config": {"max_value": 4, "axis": 0}}, "bias_constraint": {"class_name": "MaxNorm", "config": {"max_value": 4, "axis": 0}}}, "name": "layer1", "inbound_nodes": [[["dropout", 0, 0, {}]]]}], "input_layers": [["input_1", 0, 0]], "output_layers": [["layer1", 0, 0]]}}, {"class_name": "Dense", "config": {"name": "dense", "trainable": true, "dtype": "float32", "units": 100, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dense", "config": {"name": "dense_1", "trainable": true, "dtype": "float32", "units": 10, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}]}, "name": "sequential", "inbound_nodes": [[["input_1", 0, 0, {}]], [["input_2", 0, 0, {}]]]}, {"class_name": "Lambda", "config": {"name": "lambda", "trainable": true, "dtype": "float32", "function": {"class_name": "__tuple__", "items": ["4wEAAAAAAAAAAQAAAAUAAABTAAAAcxYAAAB0AKABfABkARkAfABkAhkAGAChAVMAKQNO6QAAAADp\nAQAAACkC2gFL2gNhYnMpAdoHdGVuc29yc6kAcgYAAAD6K3NpYW1lc2VOZXRfbm9jb21tb25jb21w\nX2xhcmdlYmF0Y2hlc19leHQucHnaCDxsYW1iZGE+QgAAAPMAAAAA\n", null, null]}, "function_type": "lambda", "module": "__main__", "output_shape": null, "output_shape_type": "raw", "output_shape_module": null, "arguments": {}}, "name": "lambda", "inbound_nodes": [[["sequential", 1, 0, {}], ["sequential", 2, 0, {}]]]}, {"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 10]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_3"}, "name": "input_3", "inbound_nodes": []}, {"class_name": "Concatenate", "config": {"name": "concatenate", "trainable": true, "dtype": "float32", "axis": -1}, "name": "concatenate", "inbound_nodes": [[["lambda", 0, 0, {}], ["input_3", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "dense_2", "trainable": true, "dtype": "float32", "units": 1, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "dense_2", "inbound_nodes": [[["concatenate", 0, 0, {}]]]}], "input_layers": [["input_1", 0, 0], ["input_2", 0, 0], ["input_3", 0, 0]], "output_layers": [["dense_2", 0, 0]]}, "build_input_shape": [{"class_name": "TensorShape", "items": [null, 24, 24, 24, 167]}, {"class_name": "TensorShape", "items": [null, 24, 24, 24, 167]}, {"class_name": "TensorShape", "items": [null, 10]}], "is_graph_network": true, "keras_version": "2.4.0", "backend": "tensorflow", "model_config": {"class_name": "Functional", "config": {"name": "functional_3", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 24, 24, 24, 167]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_1"}, "name": "input_1", "inbound_nodes": []}, {"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 24, 24, 24, 167]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_2"}, "name": "input_2", "inbound_nodes": []}, {"class_name": "Sequential", "config": {"name": "sequential", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 24, 24, 24, 167]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "functional_1_input"}}, {"class_name": "Functional", "config": {"name": "functional_1", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 24, 24, 24, 167]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_1"}, "name": "input_1", "inbound_nodes": []}, {"class_name": "Conv3D", "config": {"name": "conv3d", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 24, 24, 24, 167]}, "dtype": "float32", "filters": 20, "kernel_size": {"class_name": "__tuple__", "items": [1, 1, 1]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d", "inbound_nodes": [[["input_1", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization", "inbound_nodes": [[["conv3d", 0, 0, {}]]]}, {"class_name": "Conv3D", "config": {"name": "conv3d_1", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 24, 24, 24, 20]}, "dtype": "float32", "filters": 20, "kernel_size": {"class_name": "__tuple__", "items": [3, 3, 3]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d_1", "inbound_nodes": [[["batch_normalization", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization_1", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization_1", "inbound_nodes": [[["conv3d_1", 0, 0, {}]]]}, {"class_name": "Conv3D", "config": {"name": "conv3d_2", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 22, 22, 22, 20]}, "dtype": "float32", "filters": 30, "kernel_size": {"class_name": "__tuple__", "items": [4, 4, 4]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d_2", "inbound_nodes": [[["batch_normalization_1", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization_2", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization_2", "inbound_nodes": [[["conv3d_2", 0, 0, {}]]]}, {"class_name": "Conv3D", "config": {"name": "conv3d_3", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 19, 19, 19, 30]}, "dtype": "float32", "filters": 20, "kernel_size": {"class_name": "__tuple__", "items": [4, 4, 4]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d_3", "inbound_nodes": [[["batch_normalization_2", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization_3", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization_3", "inbound_nodes": [[["conv3d_3", 0, 0, {}]]]}, {"class_name": "AveragePooling3D", "config": {"name": "average_pooling3d", "trainable": true, "dtype": "float32", "pool_size": {"class_name": "__tuple__", "items": [4, 4, 4]}, "padding": "valid", "strides": {"class_name": "__tuple__", "items": [4, 4, 4]}, "data_format": "channels_last"}, "name": "average_pooling3d", "inbound_nodes": [[["batch_normalization_3", 0, 0, {}]]]}, {"class_name": "Flatten", "config": {"name": "flatten", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "name": "flatten", "inbound_nodes": [[["average_pooling3d", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "dropout", "trainable": true, "dtype": "float32", "rate": 0.4, "noise_shape": null, "seed": null}, "name": "dropout", "inbound_nodes": [[["flatten", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "layer1", "trainable": true, "dtype": "float32", "units": 200, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": {"class_name": "MaxNorm", "config": {"max_value": 4, "axis": 0}}, "bias_constraint": {"class_name": "MaxNorm", "config": {"max_value": 4, "axis": 0}}}, "name": "layer1", "inbound_nodes": [[["dropout", 0, 0, {}]]]}], "input_layers": [["input_1", 0, 0]], "output_layers": [["layer1", 0, 0]]}}, {"class_name": "Dense", "config": {"name": "dense", "trainable": true, "dtype": "float32", "units": 100, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dense", "config": {"name": "dense_1", "trainable": true, "dtype": "float32", "units": 10, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}]}, "name": "sequential", "inbound_nodes": [[["input_1", 0, 0, {}]], [["input_2", 0, 0, {}]]]}, {"class_name": "Lambda", "config": {"name": "lambda", "trainable": true, "dtype": "float32", "function": {"class_name": "__tuple__", "items": ["4wEAAAAAAAAAAQAAAAUAAABTAAAAcxYAAAB0AKABfABkARkAfABkAhkAGAChAVMAKQNO6QAAAADp\nAQAAACkC2gFL2gNhYnMpAdoHdGVuc29yc6kAcgYAAAD6K3NpYW1lc2VOZXRfbm9jb21tb25jb21w\nX2xhcmdlYmF0Y2hlc19leHQucHnaCDxsYW1iZGE+QgAAAPMAAAAA\n", null, null]}, "function_type": "lambda", "module": "__main__", "output_shape": null, "output_shape_type": "raw", "output_shape_module": null, "arguments": {}}, "name": "lambda", "inbound_nodes": [[["sequential", 1, 0, {}], ["sequential", 2, 0, {}]]]}, {"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 10]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_3"}, "name": "input_3", "inbound_nodes": []}, {"class_name": "Concatenate", "config": {"name": "concatenate", "trainable": true, "dtype": "float32", "axis": -1}, "name": "concatenate", "inbound_nodes": [[["lambda", 0, 0, {}], ["input_3", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "dense_2", "trainable": true, "dtype": "float32", "units": 1, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "dense_2", "inbound_nodes": [[["concatenate", 0, 0, {}]]]}], "input_layers": [["input_1", 0, 0], ["input_2", 0, 0], ["input_3", 0, 0]], "output_layers": [["dense_2", 0, 0]]}}, "training_config": {"loss": "mean_squared_error", "metrics": null, "weighted_metrics": null, "loss_weights": null, "optimizer_config": {"class_name": "Adam", "config": {"name": "Adam", "learning_rate": 0.0010000000474974513, "decay": 0.0, "beta_1": 0.8999999761581421, "beta_2": 0.9990000128746033, "epsilon": 1e-07, "amsgrad": false}}}}
?"?
_tf_keras_input_layer?{"class_name": "InputLayer", "name": "input_1", "dtype": "float32", "sparse": false, "ragged": false, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 24, 24, 24, 167]}, "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 24, 24, 24, 167]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_1"}}
?"?
_tf_keras_input_layer?{"class_name": "InputLayer", "name": "input_2", "dtype": "float32", "sparse": false, "ragged": false, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 24, 24, 24, 167]}, "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 24, 24, 24, 167]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_2"}}
??
layer_with_weights-0
layer-0
layer_with_weights-1
layer-1
layer_with_weights-2
layer-2
trainable_variables
regularization_losses
	variables
	keras_api
+?&call_and_return_all_conditional_losses
?__call__"??
_tf_keras_sequential??{"class_name": "Sequential", "name": "sequential", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "must_restore_from_config": false, "config": {"name": "sequential", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 24, 24, 24, 167]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "functional_1_input"}}, {"class_name": "Functional", "config": {"name": "functional_1", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 24, 24, 24, 167]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_1"}, "name": "input_1", "inbound_nodes": []}, {"class_name": "Conv3D", "config": {"name": "conv3d", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 24, 24, 24, 167]}, "dtype": "float32", "filters": 20, "kernel_size": {"class_name": "__tuple__", "items": [1, 1, 1]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d", "inbound_nodes": [[["input_1", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization", "inbound_nodes": [[["conv3d", 0, 0, {}]]]}, {"class_name": "Conv3D", "config": {"name": "conv3d_1", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 24, 24, 24, 20]}, "dtype": "float32", "filters": 20, "kernel_size": {"class_name": "__tuple__", "items": [3, 3, 3]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d_1", "inbound_nodes": [[["batch_normalization", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization_1", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization_1", "inbound_nodes": [[["conv3d_1", 0, 0, {}]]]}, {"class_name": "Conv3D", "config": {"name": "conv3d_2", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 22, 22, 22, 20]}, "dtype": "float32", "filters": 30, "kernel_size": {"class_name": "__tuple__", "items": [4, 4, 4]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d_2", "inbound_nodes": [[["batch_normalization_1", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization_2", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization_2", "inbound_nodes": [[["conv3d_2", 0, 0, {}]]]}, {"class_name": "Conv3D", "config": {"name": "conv3d_3", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 19, 19, 19, 30]}, "dtype": "float32", "filters": 20, "kernel_size": {"class_name": "__tuple__", "items": [4, 4, 4]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d_3", "inbound_nodes": [[["batch_normalization_2", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization_3", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization_3", "inbound_nodes": [[["conv3d_3", 0, 0, {}]]]}, {"class_name": "AveragePooling3D", "config": {"name": "average_pooling3d", "trainable": true, "dtype": "float32", "pool_size": {"class_name": "__tuple__", "items": [4, 4, 4]}, "padding": "valid", "strides": {"class_name": "__tuple__", "items": [4, 4, 4]}, "data_format": "channels_last"}, "name": "average_pooling3d", "inbound_nodes": [[["batch_normalization_3", 0, 0, {}]]]}, {"class_name": "Flatten", "config": {"name": "flatten", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "name": "flatten", "inbound_nodes": [[["average_pooling3d", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "dropout", "trainable": true, "dtype": "float32", "rate": 0.4, "noise_shape": null, "seed": null}, "name": "dropout", "inbound_nodes": [[["flatten", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "layer1", "trainable": true, "dtype": "float32", "units": 200, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": {"class_name": "MaxNorm", "config": {"max_value": 4, "axis": 0}}, "bias_constraint": {"class_name": "MaxNorm", "config": {"max_value": 4, "axis": 0}}}, "name": "layer1", "inbound_nodes": [[["dropout", 0, 0, {}]]]}], "input_layers": [["input_1", 0, 0]], "output_layers": [["layer1", 0, 0]]}}, {"class_name": "Dense", "config": {"name": "dense", "trainable": true, "dtype": "float32", "units": 100, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dense", "config": {"name": "dense_1", "trainable": true, "dtype": "float32", "units": 10, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}]}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 24, 24, 24, 167]}, "is_graph_network": true, "keras_version": "2.4.0", "backend": "tensorflow", "model_config": {"class_name": "Sequential", "config": {"name": "sequential", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 24, 24, 24, 167]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "functional_1_input"}}, {"class_name": "Functional", "config": {"name": "functional_1", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 24, 24, 24, 167]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_1"}, "name": "input_1", "inbound_nodes": []}, {"class_name": "Conv3D", "config": {"name": "conv3d", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 24, 24, 24, 167]}, "dtype": "float32", "filters": 20, "kernel_size": {"class_name": "__tuple__", "items": [1, 1, 1]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d", "inbound_nodes": [[["input_1", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization", "inbound_nodes": [[["conv3d", 0, 0, {}]]]}, {"class_name": "Conv3D", "config": {"name": "conv3d_1", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 24, 24, 24, 20]}, "dtype": "float32", "filters": 20, "kernel_size": {"class_name": "__tuple__", "items": [3, 3, 3]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d_1", "inbound_nodes": [[["batch_normalization", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization_1", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization_1", "inbound_nodes": [[["conv3d_1", 0, 0, {}]]]}, {"class_name": "Conv3D", "config": {"name": "conv3d_2", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 22, 22, 22, 20]}, "dtype": "float32", "filters": 30, "kernel_size": {"class_name": "__tuple__", "items": [4, 4, 4]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d_2", "inbound_nodes": [[["batch_normalization_1", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization_2", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization_2", "inbound_nodes": [[["conv3d_2", 0, 0, {}]]]}, {"class_name": "Conv3D", "config": {"name": "conv3d_3", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 19, 19, 19, 30]}, "dtype": "float32", "filters": 20, "kernel_size": {"class_name": "__tuple__", "items": [4, 4, 4]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d_3", "inbound_nodes": [[["batch_normalization_2", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization_3", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization_3", "inbound_nodes": [[["conv3d_3", 0, 0, {}]]]}, {"class_name": "AveragePooling3D", "config": {"name": "average_pooling3d", "trainable": true, "dtype": "float32", "pool_size": {"class_name": "__tuple__", "items": [4, 4, 4]}, "padding": "valid", "strides": {"class_name": "__tuple__", "items": [4, 4, 4]}, "data_format": "channels_last"}, "name": "average_pooling3d", "inbound_nodes": [[["batch_normalization_3", 0, 0, {}]]]}, {"class_name": "Flatten", "config": {"name": "flatten", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "name": "flatten", "inbound_nodes": [[["average_pooling3d", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "dropout", "trainable": true, "dtype": "float32", "rate": 0.4, "noise_shape": null, "seed": null}, "name": "dropout", "inbound_nodes": [[["flatten", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "layer1", "trainable": true, "dtype": "float32", "units": 200, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": {"class_name": "MaxNorm", "config": {"max_value": 4, "axis": 0}}, "bias_constraint": {"class_name": "MaxNorm", "config": {"max_value": 4, "axis": 0}}}, "name": "layer1", "inbound_nodes": [[["dropout", 0, 0, {}]]]}], "input_layers": [["input_1", 0, 0]], "output_layers": [["layer1", 0, 0]]}}, {"class_name": "Dense", "config": {"name": "dense", "trainable": true, "dtype": "float32", "units": 100, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}, {"class_name": "Dense", "config": {"name": "dense_1", "trainable": true, "dtype": "float32", "units": 10, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}]}}}
?
trainable_variables
regularization_losses
	variables
	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Lambda", "name": "lambda", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "lambda", "trainable": true, "dtype": "float32", "function": {"class_name": "__tuple__", "items": ["4wEAAAAAAAAAAQAAAAUAAABTAAAAcxYAAAB0AKABfABkARkAfABkAhkAGAChAVMAKQNO6QAAAADp\nAQAAACkC2gFL2gNhYnMpAdoHdGVuc29yc6kAcgYAAAD6K3NpYW1lc2VOZXRfbm9jb21tb25jb21w\nX2xhcmdlYmF0Y2hlc19leHQucHnaCDxsYW1iZGE+QgAAAPMAAAAA\n", null, null]}, "function_type": "lambda", "module": "__main__", "output_shape": null, "output_shape_type": "raw", "output_shape_module": null, "arguments": {}}}
?"?
_tf_keras_input_layer?{"class_name": "InputLayer", "name": "input_3", "dtype": "float32", "sparse": false, "ragged": false, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 10]}, "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 10]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_3"}}
?
trainable_variables
regularization_losses
	variables
	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Concatenate", "name": "concatenate", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "concatenate", "trainable": true, "dtype": "float32", "axis": -1}, "build_input_shape": [{"class_name": "TensorShape", "items": [null, 10]}, {"class_name": "TensorShape", "items": [null, 10]}]}
?

kernel
bias
trainable_variables
 regularization_losses
!	variables
"	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Dense", "name": "dense_2", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "dense_2", "trainable": true, "dtype": "float32", "units": 1, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 20}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 20]}}
?
#iter

$beta_1

%beta_2
	&decay
'learning_ratem?m?(m?)m?*m?+m?,m?-m?.m?/m?0m?1m?2m?3m?4m?5m?6m?7m?8m?9m?:m?;m?<m?=m?v?v?(v?)v?*v?+v?,v?-v?.v?/v?0v?1v?2v?3v?4v?5v?6v?7v?8v?9v?:v?;v?<v?=v?"
	optimizer
?
(0
)1
*2
+3
,4
-5
.6
/7
08
19
210
311
412
513
614
715
816
917
:18
;19
<20
=21
22
23"
trackable_list_wrapper
 "
trackable_list_wrapper
?
(0
)1
*2
+3
>4
?5
,6
-7
.8
/9
@10
A11
012
113
214
315
B16
C17
418
519
620
721
D22
E23
824
925
:26
;27
<28
=29
30
31"
trackable_list_wrapper
?
	trainable_variables

regularization_losses
	variables

Flayers
Glayer_metrics
Hmetrics
Inon_trainable_variables
Jlayer_regularization_losses
?__call__
?_default_save_signature
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
-
?serving_default"
signature_map
?
Klayer-0
Llayer_with_weights-0
Llayer-1
Mlayer_with_weights-1
Mlayer-2
Nlayer_with_weights-2
Nlayer-3
Olayer_with_weights-3
Olayer-4
Player_with_weights-4
Player-5
Qlayer_with_weights-5
Qlayer-6
Rlayer_with_weights-6
Rlayer-7
Slayer_with_weights-7
Slayer-8
Tlayer-9
Ulayer-10
Vlayer-11
Wlayer_with_weights-8
Wlayer-12
Xtrainable_variables
Yregularization_losses
Z	variables
[	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?z
_tf_keras_network?z{"class_name": "Functional", "name": "functional_1", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "must_restore_from_config": false, "config": {"name": "functional_1", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 24, 24, 24, 167]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_1"}, "name": "input_1", "inbound_nodes": []}, {"class_name": "Conv3D", "config": {"name": "conv3d", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 24, 24, 24, 167]}, "dtype": "float32", "filters": 20, "kernel_size": {"class_name": "__tuple__", "items": [1, 1, 1]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d", "inbound_nodes": [[["input_1", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization", "inbound_nodes": [[["conv3d", 0, 0, {}]]]}, {"class_name": "Conv3D", "config": {"name": "conv3d_1", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 24, 24, 24, 20]}, "dtype": "float32", "filters": 20, "kernel_size": {"class_name": "__tuple__", "items": [3, 3, 3]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d_1", "inbound_nodes": [[["batch_normalization", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization_1", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization_1", "inbound_nodes": [[["conv3d_1", 0, 0, {}]]]}, {"class_name": "Conv3D", "config": {"name": "conv3d_2", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 22, 22, 22, 20]}, "dtype": "float32", "filters": 30, "kernel_size": {"class_name": "__tuple__", "items": [4, 4, 4]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d_2", "inbound_nodes": [[["batch_normalization_1", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization_2", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization_2", "inbound_nodes": [[["conv3d_2", 0, 0, {}]]]}, {"class_name": "Conv3D", "config": {"name": "conv3d_3", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 19, 19, 19, 30]}, "dtype": "float32", "filters": 20, "kernel_size": {"class_name": "__tuple__", "items": [4, 4, 4]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d_3", "inbound_nodes": [[["batch_normalization_2", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization_3", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization_3", "inbound_nodes": [[["conv3d_3", 0, 0, {}]]]}, {"class_name": "AveragePooling3D", "config": {"name": "average_pooling3d", "trainable": true, "dtype": "float32", "pool_size": {"class_name": "__tuple__", "items": [4, 4, 4]}, "padding": "valid", "strides": {"class_name": "__tuple__", "items": [4, 4, 4]}, "data_format": "channels_last"}, "name": "average_pooling3d", "inbound_nodes": [[["batch_normalization_3", 0, 0, {}]]]}, {"class_name": "Flatten", "config": {"name": "flatten", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "name": "flatten", "inbound_nodes": [[["average_pooling3d", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "dropout", "trainable": true, "dtype": "float32", "rate": 0.4, "noise_shape": null, "seed": null}, "name": "dropout", "inbound_nodes": [[["flatten", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "layer1", "trainable": true, "dtype": "float32", "units": 200, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": {"class_name": "MaxNorm", "config": {"max_value": 4, "axis": 0}}, "bias_constraint": {"class_name": "MaxNorm", "config": {"max_value": 4, "axis": 0}}}, "name": "layer1", "inbound_nodes": [[["dropout", 0, 0, {}]]]}], "input_layers": [["input_1", 0, 0]], "output_layers": [["layer1", 0, 0]]}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 24, 24, 24, 167]}, "is_graph_network": true, "keras_version": "2.4.0", "backend": "tensorflow", "model_config": {"class_name": "Functional", "config": {"name": "functional_1", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 24, 24, 24, 167]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_1"}, "name": "input_1", "inbound_nodes": []}, {"class_name": "Conv3D", "config": {"name": "conv3d", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 24, 24, 24, 167]}, "dtype": "float32", "filters": 20, "kernel_size": {"class_name": "__tuple__", "items": [1, 1, 1]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d", "inbound_nodes": [[["input_1", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization", "inbound_nodes": [[["conv3d", 0, 0, {}]]]}, {"class_name": "Conv3D", "config": {"name": "conv3d_1", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 24, 24, 24, 20]}, "dtype": "float32", "filters": 20, "kernel_size": {"class_name": "__tuple__", "items": [3, 3, 3]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d_1", "inbound_nodes": [[["batch_normalization", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization_1", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization_1", "inbound_nodes": [[["conv3d_1", 0, 0, {}]]]}, {"class_name": "Conv3D", "config": {"name": "conv3d_2", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 22, 22, 22, 20]}, "dtype": "float32", "filters": 30, "kernel_size": {"class_name": "__tuple__", "items": [4, 4, 4]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d_2", "inbound_nodes": [[["batch_normalization_1", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization_2", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization_2", "inbound_nodes": [[["conv3d_2", 0, 0, {}]]]}, {"class_name": "Conv3D", "config": {"name": "conv3d_3", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 19, 19, 19, 30]}, "dtype": "float32", "filters": 20, "kernel_size": {"class_name": "__tuple__", "items": [4, 4, 4]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "name": "conv3d_3", "inbound_nodes": [[["batch_normalization_2", 0, 0, {}]]]}, {"class_name": "BatchNormalization", "config": {"name": "batch_normalization_3", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "name": "batch_normalization_3", "inbound_nodes": [[["conv3d_3", 0, 0, {}]]]}, {"class_name": "AveragePooling3D", "config": {"name": "average_pooling3d", "trainable": true, "dtype": "float32", "pool_size": {"class_name": "__tuple__", "items": [4, 4, 4]}, "padding": "valid", "strides": {"class_name": "__tuple__", "items": [4, 4, 4]}, "data_format": "channels_last"}, "name": "average_pooling3d", "inbound_nodes": [[["batch_normalization_3", 0, 0, {}]]]}, {"class_name": "Flatten", "config": {"name": "flatten", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "name": "flatten", "inbound_nodes": [[["average_pooling3d", 0, 0, {}]]]}, {"class_name": "Dropout", "config": {"name": "dropout", "trainable": true, "dtype": "float32", "rate": 0.4, "noise_shape": null, "seed": null}, "name": "dropout", "inbound_nodes": [[["flatten", 0, 0, {}]]]}, {"class_name": "Dense", "config": {"name": "layer1", "trainable": true, "dtype": "float32", "units": 200, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": {"class_name": "MaxNorm", "config": {"max_value": 4, "axis": 0}}, "bias_constraint": {"class_name": "MaxNorm", "config": {"max_value": 4, "axis": 0}}}, "name": "layer1", "inbound_nodes": [[["dropout", 0, 0, {}]]]}], "input_layers": [["input_1", 0, 0]], "output_layers": [["layer1", 0, 0]]}}}
?

:kernel
;bias
\trainable_variables
]regularization_losses
^	variables
_	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Dense", "name": "dense", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "dense", "trainable": true, "dtype": "float32", "units": 100, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 200}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 200]}}
?

<kernel
=bias
`trainable_variables
aregularization_losses
b	variables
c	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Dense", "name": "dense_1", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "dense_1", "trainable": true, "dtype": "float32", "units": 10, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": {"class_name": "L2", "config": {"l2": 0.0010000000474974513}}, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 100}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 100]}}
?
(0
)1
*2
+3
,4
-5
.6
/7
08
19
210
311
412
513
614
715
816
917
:18
;19
<20
=21"
trackable_list_wrapper
0
?0
?1"
trackable_list_wrapper
?
(0
)1
*2
+3
>4
?5
,6
-7
.8
/9
@10
A11
012
113
214
315
B16
C17
418
519
620
721
D22
E23
824
925
:26
;27
<28
=29"
trackable_list_wrapper
?
trainable_variables
regularization_losses
	variables

dlayers
elayer_metrics
fmetrics
gnon_trainable_variables
hlayer_regularization_losses
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
trainable_variables
regularization_losses
	variables
ilayer_metrics

jlayers
kmetrics
lnon_trainable_variables
mlayer_regularization_losses
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
trainable_variables
regularization_losses
	variables
nlayer_metrics

olayers
pmetrics
qnon_trainable_variables
rlayer_regularization_losses
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 :2dense_2/kernel
:2dense_2/bias
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
?
trainable_variables
 regularization_losses
!	variables
slayer_metrics

tlayers
umetrics
vnon_trainable_variables
wlayer_regularization_losses
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
:	 (2	Adam/iter
: (2Adam/beta_1
: (2Adam/beta_2
: (2
Adam/decay
: (2Adam/learning_rate
,:*?2conv3d/kernel
:2conv3d/bias
':%2batch_normalization/gamma
&:$2batch_normalization/beta
-:+2conv3d_1/kernel
:2conv3d_1/bias
):'2batch_normalization_1/gamma
(:&2batch_normalization_1/beta
-:+2conv3d_2/kernel
:2conv3d_2/bias
):'2batch_normalization_2/gamma
(:&2batch_normalization_2/beta
-:+2conv3d_3/kernel
:2conv3d_3/bias
):'2batch_normalization_3/gamma
(:&2batch_normalization_3/beta
!:
?
?2layer1/kernel
:?2layer1/bias
:	?d2dense/kernel
:d2
dense/bias
 :d
2dense_1/kernel
:
2dense_1/bias
/:- (2batch_normalization/moving_mean
3:1 (2#batch_normalization/moving_variance
1:/ (2!batch_normalization_1/moving_mean
5:3 (2%batch_normalization_1/moving_variance
1:/ (2!batch_normalization_2/moving_mean
5:3 (2%batch_normalization_2/moving_variance
1:/ (2!batch_normalization_3/moving_mean
5:3 (2%batch_normalization_3/moving_variance
Q
0
1
2
3
4
5
6"
trackable_list_wrapper
 "
trackable_dict_wrapper
'
x0"
trackable_list_wrapper
X
>0
?1
@2
A3
B4
C5
D6
E7"
trackable_list_wrapper
 "
trackable_list_wrapper
?
#y_self_saveable_object_factories"?
_tf_keras_input_layer?{"class_name": "InputLayer", "name": "input_1", "dtype": "float32", "sparse": false, "ragged": false, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 24, 24, 24, 167]}, "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 24, 24, 24, 167]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "input_1"}}
?

(kernel
)bias
#z_self_saveable_object_factories
{trainable_variables
|regularization_losses
}	variables
~	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?	
_tf_keras_layer?	{"class_name": "Conv3D", "name": "conv3d", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 24, 24, 24, 167]}, "stateful": false, "must_restore_from_config": false, "config": {"name": "conv3d", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 24, 24, 24, 167]}, "dtype": "float32", "filters": 20, "kernel_size": {"class_name": "__tuple__", "items": [1, 1, 1]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "linear", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 5, "axes": {"-1": 167}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 24, 24, 24, 167]}}
?	
axis
	*gamma
+beta
>moving_mean
?moving_variance
$?_self_saveable_object_factories
?trainable_variables
?regularization_losses
?	variables
?	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "BatchNormalization", "name": "batch_normalization", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "batch_normalization", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 5, "max_ndim": null, "min_ndim": null, "axes": {"4": 20}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 24, 24, 24, 20]}}
?

,kernel
-bias
$?_self_saveable_object_factories
?trainable_variables
?regularization_losses
?	variables
?	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?	
_tf_keras_layer?	{"class_name": "Conv3D", "name": "conv3d_1", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 24, 24, 24, 20]}, "stateful": false, "must_restore_from_config": false, "config": {"name": "conv3d_1", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 24, 24, 24, 20]}, "dtype": "float32", "filters": 20, "kernel_size": {"class_name": "__tuple__", "items": [3, 3, 3]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 5, "axes": {"-1": 20}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 24, 24, 24, 20]}}
?	
	?axis
	.gamma
/beta
@moving_mean
Amoving_variance
$?_self_saveable_object_factories
?trainable_variables
?regularization_losses
?	variables
?	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "BatchNormalization", "name": "batch_normalization_1", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "batch_normalization_1", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 5, "max_ndim": null, "min_ndim": null, "axes": {"4": 20}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 22, 22, 22, 20]}}
?

0kernel
1bias
$?_self_saveable_object_factories
?trainable_variables
?regularization_losses
?	variables
?	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?	
_tf_keras_layer?	{"class_name": "Conv3D", "name": "conv3d_2", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 22, 22, 22, 20]}, "stateful": false, "must_restore_from_config": false, "config": {"name": "conv3d_2", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 22, 22, 22, 20]}, "dtype": "float32", "filters": 30, "kernel_size": {"class_name": "__tuple__", "items": [4, 4, 4]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 5, "axes": {"-1": 20}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 22, 22, 22, 20]}}
?	
	?axis
	2gamma
3beta
Bmoving_mean
Cmoving_variance
$?_self_saveable_object_factories
?trainable_variables
?regularization_losses
?	variables
?	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "BatchNormalization", "name": "batch_normalization_2", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "batch_normalization_2", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 5, "max_ndim": null, "min_ndim": null, "axes": {"4": 30}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 19, 19, 19, 30]}}
?

4kernel
5bias
$?_self_saveable_object_factories
?trainable_variables
?regularization_losses
?	variables
?	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?	
_tf_keras_layer?	{"class_name": "Conv3D", "name": "conv3d_3", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 19, 19, 19, 30]}, "stateful": false, "must_restore_from_config": false, "config": {"name": "conv3d_3", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, null, 19, 19, 19, 30]}, "dtype": "float32", "filters": 20, "kernel_size": {"class_name": "__tuple__", "items": [4, 4, 4]}, "strides": {"class_name": "__tuple__", "items": [1, 1, 1]}, "padding": "valid", "data_format": "channels_last", "dilation_rate": {"class_name": "__tuple__", "items": [1, 1, 1]}, "groups": 1, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "HeUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 5, "axes": {"-1": 30}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 19, 19, 19, 30]}}
?	
	?axis
	6gamma
7beta
Dmoving_mean
Emoving_variance
$?_self_saveable_object_factories
?trainable_variables
?regularization_losses
?	variables
?	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "BatchNormalization", "name": "batch_normalization_3", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "batch_normalization_3", "trainable": true, "dtype": "float32", "axis": [4], "momentum": 0.99, "epsilon": 0.001, "center": true, "scale": true, "beta_initializer": {"class_name": "Zeros", "config": {}}, "gamma_initializer": {"class_name": "Ones", "config": {}}, "moving_mean_initializer": {"class_name": "Zeros", "config": {}}, "moving_variance_initializer": {"class_name": "Ones", "config": {}}, "beta_regularizer": null, "gamma_regularizer": null, "beta_constraint": null, "gamma_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 5, "max_ndim": null, "min_ndim": null, "axes": {"4": 20}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 16, 16, 16, 20]}}
?
$?_self_saveable_object_factories
?trainable_variables
?regularization_losses
?	variables
?	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "AveragePooling3D", "name": "average_pooling3d", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "average_pooling3d", "trainable": true, "dtype": "float32", "pool_size": {"class_name": "__tuple__", "items": [4, 4, 4]}, "padding": "valid", "strides": {"class_name": "__tuple__", "items": [4, 4, 4]}, "data_format": "channels_last"}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": 5, "max_ndim": null, "min_ndim": null, "axes": {}}}}
?
$?_self_saveable_object_factories
?trainable_variables
?regularization_losses
?	variables
?	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Flatten", "name": "flatten", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "flatten", "trainable": true, "dtype": "float32", "data_format": "channels_last"}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 1, "axes": {}}}}
?
$?_self_saveable_object_factories
?trainable_variables
?regularization_losses
?	variables
?	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Dropout", "name": "dropout", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "dropout", "trainable": true, "dtype": "float32", "rate": 0.4, "noise_shape": null, "seed": null}}
?	

8kernel
9bias
$?_self_saveable_object_factories
?trainable_variables
?regularization_losses
?	variables
?	keras_api
+?&call_and_return_all_conditional_losses
?__call__"?
_tf_keras_layer?{"class_name": "Dense", "name": "layer1", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "layer1", "trainable": true, "dtype": "float32", "units": 200, "activation": "elu", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": {"class_name": "MaxNorm", "config": {"max_value": 4, "axis": 0}}, "bias_constraint": {"class_name": "MaxNorm", "config": {"max_value": 4, "axis": 0}}}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 1280}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 1280]}}
?
(0
)1
*2
+3
,4
-5
.6
/7
08
19
210
311
412
513
614
715
816
917"
trackable_list_wrapper
 "
trackable_list_wrapper
?
(0
)1
*2
+3
>4
?5
,6
-7
.8
/9
@10
A11
012
113
214
315
B16
C17
418
519
620
721
D22
E23
824
925"
trackable_list_wrapper
?
Xtrainable_variables
Yregularization_losses
Z	variables
?layers
?layer_metrics
?metrics
?non_trainable_variables
 ?layer_regularization_losses
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
.
:0
;1"
trackable_list_wrapper
(
?0"
trackable_list_wrapper
.
:0
;1"
trackable_list_wrapper
?
\trainable_variables
]regularization_losses
^	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
.
<0
=1"
trackable_list_wrapper
(
?0"
trackable_list_wrapper
.
<0
=1"
trackable_list_wrapper
?
`trainable_variables
aregularization_losses
b	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
5
0
1
2"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
X
>0
?1
@2
A3
B4
C5
D6
E7"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?

?total

?count
?	variables
?	keras_api"?
_tf_keras_metricj{"class_name": "Mean", "name": "loss", "dtype": "float32", "config": {"name": "loss", "dtype": "float32"}}
 "
trackable_dict_wrapper
 "
trackable_dict_wrapper
.
(0
)1"
trackable_list_wrapper
 "
trackable_list_wrapper
.
(0
)1"
trackable_list_wrapper
?
{trainable_variables
|regularization_losses
}	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
.
*0
+1"
trackable_list_wrapper
 "
trackable_list_wrapper
<
*0
+1
>2
?3"
trackable_list_wrapper
?
?trainable_variables
?regularization_losses
?	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_dict_wrapper
.
,0
-1"
trackable_list_wrapper
 "
trackable_list_wrapper
.
,0
-1"
trackable_list_wrapper
?
?trainable_variables
?regularization_losses
?	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
.
.0
/1"
trackable_list_wrapper
 "
trackable_list_wrapper
<
.0
/1
@2
A3"
trackable_list_wrapper
?
?trainable_variables
?regularization_losses
?	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_dict_wrapper
.
00
11"
trackable_list_wrapper
 "
trackable_list_wrapper
.
00
11"
trackable_list_wrapper
?
?trainable_variables
?regularization_losses
?	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
.
20
31"
trackable_list_wrapper
 "
trackable_list_wrapper
<
20
31
B2
C3"
trackable_list_wrapper
?
?trainable_variables
?regularization_losses
?	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_dict_wrapper
.
40
51"
trackable_list_wrapper
 "
trackable_list_wrapper
.
40
51"
trackable_list_wrapper
?
?trainable_variables
?regularization_losses
?	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
.
60
71"
trackable_list_wrapper
 "
trackable_list_wrapper
<
60
71
D2
E3"
trackable_list_wrapper
?
?trainable_variables
?regularization_losses
?	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
?trainable_variables
?regularization_losses
?	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
?trainable_variables
?regularization_losses
?	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
?trainable_variables
?regularization_losses
?	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_dict_wrapper
.
80
91"
trackable_list_wrapper
 "
trackable_list_wrapper
.
80
91"
trackable_list_wrapper
?
?trainable_variables
?regularization_losses
?	variables
?layer_metrics
?layers
?metrics
?non_trainable_variables
 ?layer_regularization_losses
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
~
K0
L1
M2
N3
O4
P5
Q6
R7
S8
T9
U10
V11
W12"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
X
>0
?1
@2
A3
B4
C5
D6
E7"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
(
?0"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
(
?0"
trackable_list_wrapper
:  (2total
:  (2count
0
?0
?1"
trackable_list_wrapper
.
?	variables"
_generic_user_object
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
.
>0
?1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
.
@0
A1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
.
B0
C1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
.
D0
E1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
%:#2Adam/dense_2/kernel/m
:2Adam/dense_2/bias/m
1:/?2Adam/conv3d/kernel/m
:2Adam/conv3d/bias/m
,:*2 Adam/batch_normalization/gamma/m
+:)2Adam/batch_normalization/beta/m
2:02Adam/conv3d_1/kernel/m
 :2Adam/conv3d_1/bias/m
.:,2"Adam/batch_normalization_1/gamma/m
-:+2!Adam/batch_normalization_1/beta/m
2:02Adam/conv3d_2/kernel/m
 :2Adam/conv3d_2/bias/m
.:,2"Adam/batch_normalization_2/gamma/m
-:+2!Adam/batch_normalization_2/beta/m
2:02Adam/conv3d_3/kernel/m
 :2Adam/conv3d_3/bias/m
.:,2"Adam/batch_normalization_3/gamma/m
-:+2!Adam/batch_normalization_3/beta/m
&:$
?
?2Adam/layer1/kernel/m
:?2Adam/layer1/bias/m
$:"	?d2Adam/dense/kernel/m
:d2Adam/dense/bias/m
%:#d
2Adam/dense_1/kernel/m
:
2Adam/dense_1/bias/m
%:#2Adam/dense_2/kernel/v
:2Adam/dense_2/bias/v
1:/?2Adam/conv3d/kernel/v
:2Adam/conv3d/bias/v
,:*2 Adam/batch_normalization/gamma/v
+:)2Adam/batch_normalization/beta/v
2:02Adam/conv3d_1/kernel/v
 :2Adam/conv3d_1/bias/v
.:,2"Adam/batch_normalization_1/gamma/v
-:+2!Adam/batch_normalization_1/beta/v
2:02Adam/conv3d_2/kernel/v
 :2Adam/conv3d_2/bias/v
.:,2"Adam/batch_normalization_2/gamma/v
-:+2!Adam/batch_normalization_2/beta/v
2:02Adam/conv3d_3/kernel/v
 :2Adam/conv3d_3/bias/v
.:,2"Adam/batch_normalization_3/gamma/v
-:+2!Adam/batch_normalization_3/beta/v
&:$
?
?2Adam/layer1/kernel/v
:?2Adam/layer1/bias/v
$:"	?d2Adam/dense/kernel/v
:d2Adam/dense/bias/v
%:#d
2Adam/dense_1/kernel/v
:
2Adam/dense_1/bias/v
?2?
I__inference_functional_3_layer_call_and_return_conditional_losses_1224678
I__inference_functional_3_layer_call_and_return_conditional_losses_1225750
I__inference_functional_3_layer_call_and_return_conditional_losses_1225520
I__inference_functional_3_layer_call_and_return_conditional_losses_1224561?
???
FullArgSpec1
args)?&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults?
p 

 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
.__inference_functional_3_layer_call_fn_1225055
.__inference_functional_3_layer_call_fn_1225821
.__inference_functional_3_layer_call_fn_1225892
.__inference_functional_3_layer_call_fn_1224867?
???
FullArgSpec1
args)?&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults?
p 

 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
"__inference__wrapped_model_1222192?
???
FullArgSpec
args? 
varargsjargs
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *???
???
.?+
input_1??????????
.?+
input_2??????????
!?
input_3?????????

?2?
G__inference_sequential_layer_call_and_return_conditional_losses_1226239
G__inference_sequential_layer_call_and_return_conditional_losses_1226107
G__inference_sequential_layer_call_and_return_conditional_losses_1223881
G__inference_sequential_layer_call_and_return_conditional_losses_1223960?
???
FullArgSpec1
args)?&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults?
p 

 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
,__inference_sequential_layer_call_fn_1226369
,__inference_sequential_layer_call_fn_1226304
,__inference_sequential_layer_call_fn_1224105
,__inference_sequential_layer_call_fn_1224249?
???
FullArgSpec1
args)?&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults?
p 

 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
C__inference_lambda_layer_call_and_return_conditional_losses_1226383
C__inference_lambda_layer_call_and_return_conditional_losses_1226376?
???
FullArgSpec1
args)?&
jself
jinputs
jmask

jtraining
varargs
 
varkw
 
defaults?

 
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
(__inference_lambda_layer_call_fn_1226389
(__inference_lambda_layer_call_fn_1226395?
???
FullArgSpec1
args)?&
jself
jinputs
jmask

jtraining
varargs
 
varkw
 
defaults?

 
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
H__inference_concatenate_layer_call_and_return_conditional_losses_1226402?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
-__inference_concatenate_layer_call_fn_1226408?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
D__inference_dense_2_layer_call_and_return_conditional_losses_1226418?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
)__inference_dense_2_layer_call_fn_1226427?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
DBB
%__inference_signature_wrapper_1225148input_1input_2input_3
?2?
I__inference_functional_1_layer_call_and_return_conditional_losses_1226710
I__inference_functional_1_layer_call_and_return_conditional_losses_1223380
I__inference_functional_1_layer_call_and_return_conditional_losses_1226604
I__inference_functional_1_layer_call_and_return_conditional_losses_1223312?
???
FullArgSpec1
args)?&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults?
p 

 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
.__inference_functional_1_layer_call_fn_1226767
.__inference_functional_1_layer_call_fn_1223506
.__inference_functional_1_layer_call_fn_1226824
.__inference_functional_1_layer_call_fn_1223631?
???
FullArgSpec1
args)?&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults?
p 

 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
B__inference_dense_layer_call_and_return_conditional_losses_1226847?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
'__inference_dense_layer_call_fn_1226856?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
D__inference_dense_1_layer_call_and_return_conditional_losses_1226879?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
)__inference_dense_1_layer_call_fn_1226888?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
__inference_loss_fn_0_1226899?
???
FullArgSpec
args? 
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *? 
?2?
__inference_loss_fn_1_1226910?
???
FullArgSpec
args? 
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *? 
?2?
C__inference_conv3d_layer_call_and_return_conditional_losses_1226920?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
(__inference_conv3d_layer_call_fn_1226929?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
P__inference_batch_normalization_layer_call_and_return_conditional_losses_1226985
P__inference_batch_normalization_layer_call_and_return_conditional_losses_1227047
P__inference_batch_normalization_layer_call_and_return_conditional_losses_1227067
P__inference_batch_normalization_layer_call_and_return_conditional_losses_1226965?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
5__inference_batch_normalization_layer_call_fn_1227093
5__inference_batch_normalization_layer_call_fn_1226998
5__inference_batch_normalization_layer_call_fn_1227011
5__inference_batch_normalization_layer_call_fn_1227080?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
E__inference_conv3d_1_layer_call_and_return_conditional_losses_1227104?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
*__inference_conv3d_1_layer_call_fn_1227113?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_1227149
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_1227169
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_1227251
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_1227231?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
7__inference_batch_normalization_1_layer_call_fn_1227277
7__inference_batch_normalization_1_layer_call_fn_1227182
7__inference_batch_normalization_1_layer_call_fn_1227195
7__inference_batch_normalization_1_layer_call_fn_1227264?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
E__inference_conv3d_2_layer_call_and_return_conditional_losses_1227288?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
*__inference_conv3d_2_layer_call_fn_1227297?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_1227353
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_1227415
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_1227333
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_1227435?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
7__inference_batch_normalization_2_layer_call_fn_1227379
7__inference_batch_normalization_2_layer_call_fn_1227448
7__inference_batch_normalization_2_layer_call_fn_1227366
7__inference_batch_normalization_2_layer_call_fn_1227461?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
E__inference_conv3d_3_layer_call_and_return_conditional_losses_1227472?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
*__inference_conv3d_3_layer_call_fn_1227481?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_1227517
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_1227619
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_1227537
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_1227599?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
7__inference_batch_normalization_3_layer_call_fn_1227550
7__inference_batch_normalization_3_layer_call_fn_1227645
7__inference_batch_normalization_3_layer_call_fn_1227632
7__inference_batch_normalization_3_layer_call_fn_1227563?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
N__inference_average_pooling3d_layer_call_and_return_conditional_losses_1222758?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *M?J
H?EA?????????????????????????????????????????????
?2?
3__inference_average_pooling3d_layer_call_fn_1222764?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *M?J
H?EA?????????????????????????????????????????????
?2?
D__inference_flatten_layer_call_and_return_conditional_losses_1227651?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
)__inference_flatten_layer_call_fn_1227656?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
D__inference_dropout_layer_call_and_return_conditional_losses_1227668
D__inference_dropout_layer_call_and_return_conditional_losses_1227673?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
)__inference_dropout_layer_call_fn_1227678
)__inference_dropout_layer_call_fn_1227683?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
C__inference_layer1_layer_call_and_return_conditional_losses_1227694?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
(__inference_layer1_layer_call_fn_1227703?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 ?
"__inference__wrapped_model_1222192? ()?*>+,-A.@/01C2B345E6D789:;<=???
???
???
.?+
input_1??????????
.?+
input_2??????????
!?
input_3?????????

? "1?.
,
dense_2!?
dense_2??????????
N__inference_average_pooling3d_layer_call_and_return_conditional_losses_1222758?_?\
U?R
P?M
inputsA?????????????????????????????????????????????
? "U?R
K?H
0A?????????????????????????????????????????????
? ?
3__inference_average_pooling3d_layer_call_fn_1222764?_?\
U?R
P?M
inputsA?????????????????????????????????????????????
? "H?EA??????????????????????????????????????????????
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_1227149?@A./Z?W
P?M
G?D
inputs8????????????????????????????????????
p
? "L?I
B??
08????????????????????????????????????
? ?
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_1227169?A.@/Z?W
P?M
G?D
inputs8????????????????????????????????????
p 
? "L?I
B??
08????????????????????????????????????
? ?
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_1227231z@A./??<
5?2
,?)
inputs?????????
p
? "1?.
'?$
0?????????
? ?
R__inference_batch_normalization_1_layer_call_and_return_conditional_losses_1227251zA.@/??<
5?2
,?)
inputs?????????
p 
? "1?.
'?$
0?????????
? ?
7__inference_batch_normalization_1_layer_call_fn_1227182?@A./Z?W
P?M
G?D
inputs8????????????????????????????????????
p
? "??<8?????????????????????????????????????
7__inference_batch_normalization_1_layer_call_fn_1227195?A.@/Z?W
P?M
G?D
inputs8????????????????????????????????????
p 
? "??<8?????????????????????????????????????
7__inference_batch_normalization_1_layer_call_fn_1227264m@A./??<
5?2
,?)
inputs?????????
p
? "$?!??????????
7__inference_batch_normalization_1_layer_call_fn_1227277mA.@/??<
5?2
,?)
inputs?????????
p 
? "$?!??????????
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_1227333?BC23Z?W
P?M
G?D
inputs8????????????????????????????????????
p
? "L?I
B??
08????????????????????????????????????
? ?
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_1227353?C2B3Z?W
P?M
G?D
inputs8????????????????????????????????????
p 
? "L?I
B??
08????????????????????????????????????
? ?
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_1227415zBC23??<
5?2
,?)
inputs?????????
p
? "1?.
'?$
0?????????
? ?
R__inference_batch_normalization_2_layer_call_and_return_conditional_losses_1227435zC2B3??<
5?2
,?)
inputs?????????
p 
? "1?.
'?$
0?????????
? ?
7__inference_batch_normalization_2_layer_call_fn_1227366?BC23Z?W
P?M
G?D
inputs8????????????????????????????????????
p
? "??<8?????????????????????????????????????
7__inference_batch_normalization_2_layer_call_fn_1227379?C2B3Z?W
P?M
G?D
inputs8????????????????????????????????????
p 
? "??<8?????????????????????????????????????
7__inference_batch_normalization_2_layer_call_fn_1227448mBC23??<
5?2
,?)
inputs?????????
p
? "$?!??????????
7__inference_batch_normalization_2_layer_call_fn_1227461mC2B3??<
5?2
,?)
inputs?????????
p 
? "$?!??????????
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_1227517?DE67Z?W
P?M
G?D
inputs8????????????????????????????????????
p
? "L?I
B??
08????????????????????????????????????
? ?
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_1227537?E6D7Z?W
P?M
G?D
inputs8????????????????????????????????????
p 
? "L?I
B??
08????????????????????????????????????
? ?
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_1227599zDE67??<
5?2
,?)
inputs?????????
p
? "1?.
'?$
0?????????
? ?
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_1227619zE6D7??<
5?2
,?)
inputs?????????
p 
? "1?.
'?$
0?????????
? ?
7__inference_batch_normalization_3_layer_call_fn_1227550?DE67Z?W
P?M
G?D
inputs8????????????????????????????????????
p
? "??<8?????????????????????????????????????
7__inference_batch_normalization_3_layer_call_fn_1227563?E6D7Z?W
P?M
G?D
inputs8????????????????????????????????????
p 
? "??<8?????????????????????????????????????
7__inference_batch_normalization_3_layer_call_fn_1227632mDE67??<
5?2
,?)
inputs?????????
p
? "$?!??????????
7__inference_batch_normalization_3_layer_call_fn_1227645mE6D7??<
5?2
,?)
inputs?????????
p 
? "$?!??????????
P__inference_batch_normalization_layer_call_and_return_conditional_losses_1226965?>?*+Z?W
P?M
G?D
inputs8????????????????????????????????????
p
? "L?I
B??
08????????????????????????????????????
? ?
P__inference_batch_normalization_layer_call_and_return_conditional_losses_1226985??*>+Z?W
P?M
G?D
inputs8????????????????????????????????????
p 
? "L?I
B??
08????????????????????????????????????
? ?
P__inference_batch_normalization_layer_call_and_return_conditional_losses_1227047z>?*+??<
5?2
,?)
inputs?????????
p
? "1?.
'?$
0?????????
? ?
P__inference_batch_normalization_layer_call_and_return_conditional_losses_1227067z?*>+??<
5?2
,?)
inputs?????????
p 
? "1?.
'?$
0?????????
? ?
5__inference_batch_normalization_layer_call_fn_1226998?>?*+Z?W
P?M
G?D
inputs8????????????????????????????????????
p
? "??<8?????????????????????????????????????
5__inference_batch_normalization_layer_call_fn_1227011??*>+Z?W
P?M
G?D
inputs8????????????????????????????????????
p 
? "??<8?????????????????????????????????????
5__inference_batch_normalization_layer_call_fn_1227080m>?*+??<
5?2
,?)
inputs?????????
p
? "$?!??????????
5__inference_batch_normalization_layer_call_fn_1227093m?*>+??<
5?2
,?)
inputs?????????
p 
? "$?!??????????
H__inference_concatenate_layer_call_and_return_conditional_losses_1226402?Z?W
P?M
K?H
"?
inputs/0?????????

"?
inputs/1?????????

? "%?"
?
0?????????
? ?
-__inference_concatenate_layer_call_fn_1226408vZ?W
P?M
K?H
"?
inputs/0?????????

"?
inputs/1?????????

? "???????????
E__inference_conv3d_1_layer_call_and_return_conditional_losses_1227104t,-;?8
1?.
,?)
inputs?????????
? "1?.
'?$
0?????????
? ?
*__inference_conv3d_1_layer_call_fn_1227113g,-;?8
1?.
,?)
inputs?????????
? "$?!??????????
E__inference_conv3d_2_layer_call_and_return_conditional_losses_1227288t01;?8
1?.
,?)
inputs?????????
? "1?.
'?$
0?????????
? ?
*__inference_conv3d_2_layer_call_fn_1227297g01;?8
1?.
,?)
inputs?????????
? "$?!??????????
E__inference_conv3d_3_layer_call_and_return_conditional_losses_1227472t45;?8
1?.
,?)
inputs?????????
? "1?.
'?$
0?????????
? ?
*__inference_conv3d_3_layer_call_fn_1227481g45;?8
1?.
,?)
inputs?????????
? "$?!??????????
C__inference_conv3d_layer_call_and_return_conditional_losses_1226920u()<?9
2?/
-?*
inputs??????????
? "1?.
'?$
0?????????
? ?
(__inference_conv3d_layer_call_fn_1226929h()<?9
2?/
-?*
inputs??????????
? "$?!??????????
D__inference_dense_1_layer_call_and_return_conditional_losses_1226879\<=/?,
%?"
 ?
inputs?????????d
? "%?"
?
0?????????

? |
)__inference_dense_1_layer_call_fn_1226888O<=/?,
%?"
 ?
inputs?????????d
? "??????????
?
D__inference_dense_2_layer_call_and_return_conditional_losses_1226418\/?,
%?"
 ?
inputs?????????
? "%?"
?
0?????????
? |
)__inference_dense_2_layer_call_fn_1226427O/?,
%?"
 ?
inputs?????????
? "???????????
B__inference_dense_layer_call_and_return_conditional_losses_1226847]:;0?-
&?#
!?
inputs??????????
? "%?"
?
0?????????d
? {
'__inference_dense_layer_call_fn_1226856P:;0?-
&?#
!?
inputs??????????
? "??????????d?
D__inference_dropout_layer_call_and_return_conditional_losses_1227668^4?1
*?'
!?
inputs??????????

p
? "&?#
?
0??????????

? ?
D__inference_dropout_layer_call_and_return_conditional_losses_1227673^4?1
*?'
!?
inputs??????????

p 
? "&?#
?
0??????????

? ~
)__inference_dropout_layer_call_fn_1227678Q4?1
*?'
!?
inputs??????????

p
? "???????????
~
)__inference_dropout_layer_call_fn_1227683Q4?1
*?'
!?
inputs??????????

p 
? "???????????
?
D__inference_flatten_layer_call_and_return_conditional_losses_1227651e;?8
1?.
,?)
inputs?????????
? "&?#
?
0??????????

? ?
)__inference_flatten_layer_call_fn_1227656X;?8
1?.
,?)
inputs?????????
? "???????????
?
I__inference_functional_1_layer_call_and_return_conditional_losses_1223312?()>?*+,-@A./01BC2345DE6789E?B
;?8
.?+
input_1??????????
p

 
? "&?#
?
0??????????
? ?
I__inference_functional_1_layer_call_and_return_conditional_losses_1223380?()?*>+,-A.@/01C2B345E6D789E?B
;?8
.?+
input_1??????????
p 

 
? "&?#
?
0??????????
? ?
I__inference_functional_1_layer_call_and_return_conditional_losses_1226604?()>?*+,-@A./01BC2345DE6789D?A
:?7
-?*
inputs??????????
p

 
? "&?#
?
0??????????
? ?
I__inference_functional_1_layer_call_and_return_conditional_losses_1226710?()?*>+,-A.@/01C2B345E6D789D?A
:?7
-?*
inputs??????????
p 

 
? "&?#
?
0??????????
? ?
.__inference_functional_1_layer_call_fn_1223506~()>?*+,-@A./01BC2345DE6789E?B
;?8
.?+
input_1??????????
p

 
? "????????????
.__inference_functional_1_layer_call_fn_1223631~()?*>+,-A.@/01C2B345E6D789E?B
;?8
.?+
input_1??????????
p 

 
? "????????????
.__inference_functional_1_layer_call_fn_1226767}()>?*+,-@A./01BC2345DE6789D?A
:?7
-?*
inputs??????????
p

 
? "????????????
.__inference_functional_1_layer_call_fn_1226824}()?*>+,-A.@/01C2B345E6D789D?A
:?7
-?*
inputs??????????
p 

 
? "????????????
I__inference_functional_3_layer_call_and_return_conditional_losses_1224561? ()>?*+,-@A./01BC2345DE6789:;<=???
???
???
.?+
input_1??????????
.?+
input_2??????????
!?
input_3?????????

p

 
? "%?"
?
0?????????
? ?
I__inference_functional_3_layer_call_and_return_conditional_losses_1224678? ()?*>+,-A.@/01C2B345E6D789:;<=???
???
???
.?+
input_1??????????
.?+
input_2??????????
!?
input_3?????????

p 

 
? "%?"
?
0?????????
? ?
I__inference_functional_3_layer_call_and_return_conditional_losses_1225520? ()>?*+,-@A./01BC2345DE6789:;<=???
???
???
/?,
inputs/0??????????
/?,
inputs/1??????????
"?
inputs/2?????????

p

 
? "%?"
?
0?????????
? ?
I__inference_functional_3_layer_call_and_return_conditional_losses_1225750? ()?*>+,-A.@/01C2B345E6D789:;<=???
???
???
/?,
inputs/0??????????
/?,
inputs/1??????????
"?
inputs/2?????????

p 

 
? "%?"
?
0?????????
? ?
.__inference_functional_3_layer_call_fn_1224867? ()>?*+,-@A./01BC2345DE6789:;<=???
???
???
.?+
input_1??????????
.?+
input_2??????????
!?
input_3?????????

p

 
? "???????????
.__inference_functional_3_layer_call_fn_1225055? ()?*>+,-A.@/01C2B345E6D789:;<=???
???
???
.?+
input_1??????????
.?+
input_2??????????
!?
input_3?????????

p 

 
? "???????????
.__inference_functional_3_layer_call_fn_1225821? ()>?*+,-@A./01BC2345DE6789:;<=???
???
???
/?,
inputs/0??????????
/?,
inputs/1??????????
"?
inputs/2?????????

p

 
? "???????????
.__inference_functional_3_layer_call_fn_1225892? ()?*>+,-A.@/01C2B345E6D789:;<=???
???
???
/?,
inputs/0??????????
/?,
inputs/1??????????
"?
inputs/2?????????

p 

 
? "???????????
C__inference_lambda_layer_call_and_return_conditional_losses_1226376?b?_
X?U
K?H
"?
inputs/0?????????

"?
inputs/1?????????


 
p
? "%?"
?
0?????????

? ?
C__inference_lambda_layer_call_and_return_conditional_losses_1226383?b?_
X?U
K?H
"?
inputs/0?????????

"?
inputs/1?????????


 
p 
? "%?"
?
0?????????

? ?
(__inference_lambda_layer_call_fn_1226389~b?_
X?U
K?H
"?
inputs/0?????????

"?
inputs/1?????????


 
p
? "??????????
?
(__inference_lambda_layer_call_fn_1226395~b?_
X?U
K?H
"?
inputs/0?????????

"?
inputs/1?????????


 
p 
? "??????????
?
C__inference_layer1_layer_call_and_return_conditional_losses_1227694^890?-
&?#
!?
inputs??????????

? "&?#
?
0??????????
? }
(__inference_layer1_layer_call_fn_1227703Q890?-
&?#
!?
inputs??????????

? "???????????<
__inference_loss_fn_0_1226899:?

? 
? "? <
__inference_loss_fn_1_1226910<?

? 
? "? ?
G__inference_sequential_layer_call_and_return_conditional_losses_1223881?()>?*+,-@A./01BC2345DE6789:;<=P?M
F?C
9?6
functional_1_input??????????
p

 
? "%?"
?
0?????????

? ?
G__inference_sequential_layer_call_and_return_conditional_losses_1223960?()?*>+,-A.@/01C2B345E6D789:;<=P?M
F?C
9?6
functional_1_input??????????
p 

 
? "%?"
?
0?????????

? ?
G__inference_sequential_layer_call_and_return_conditional_losses_1226107?()>?*+,-@A./01BC2345DE6789:;<=D?A
:?7
-?*
inputs??????????
p

 
? "%?"
?
0?????????

? ?
G__inference_sequential_layer_call_and_return_conditional_losses_1226239?()?*>+,-A.@/01C2B345E6D789:;<=D?A
:?7
-?*
inputs??????????
p 

 
? "%?"
?
0?????????

? ?
,__inference_sequential_layer_call_fn_1224105?()>?*+,-@A./01BC2345DE6789:;<=P?M
F?C
9?6
functional_1_input??????????
p

 
? "??????????
?
,__inference_sequential_layer_call_fn_1224249?()?*>+,-A.@/01C2B345E6D789:;<=P?M
F?C
9?6
functional_1_input??????????
p 

 
? "??????????
?
,__inference_sequential_layer_call_fn_1226304?()>?*+,-@A./01BC2345DE6789:;<=D?A
:?7
-?*
inputs??????????
p

 
? "??????????
?
,__inference_sequential_layer_call_fn_1226369?()?*>+,-A.@/01C2B345E6D789:;<=D?A
:?7
-?*
inputs??????????
p 

 
? "??????????
?
%__inference_signature_wrapper_1225148? ()?*>+,-A.@/01C2B345E6D789:;<=???
? 
???
9
input_1.?+
input_1??????????
9
input_2.?+
input_2??????????
,
input_3!?
input_3?????????
"1?.
,
dense_2!?
dense_2?????????