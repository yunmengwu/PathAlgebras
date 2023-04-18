
--------------------
-- Documentation  --
--------------------
undocumented {(NewFromMethod,PathAlgebra,List),
	      (NewFromMethod,PathAlgebraQuotient,List),
	      (NewFromMethod,PAModule,List),
	      (expression, PAMatrix),
	      (numColumns, PAMatrix),
	      (numRows, PAMatrix),
	      (ideal, PathAlgebra),
	      (ideal, PathAlgebraQuotient),
	      (putInPathAlgebra,PathAlgebraQuotient,PAPath,QQ),
 	      (putInPathAlgebra,PathAlgebraQuotient,PAPath,ZZ),
	      (symbol ^, PathAlgebra, ZZ),
	      (symbol *,PAModMon,PAPath),
	      toUniform,
	      (toUniform,PAElement),
	      (toUniform,PAVector),
	      toRationalFunction,    
	      (toRationalFunction,List),
	      paHilbertSeries,
	      (paHilbertSeries,ZZ,PathAlgebra,List),	
	      paPath,
	      (paPath,PAGraph, ZZ),
	      (paPath,PAGraph,List),
	      getUnsubVar,
	      (getUnsubVar,IndexedVariable),
	      (getUnsubVar,Symbol),
	      findRecurrence,
	      (findRecurrence,List),
	      (degrees, PAModule),	      
	      freePAModule,
	      (freePAModule,PathAlgebra,ZZ),
	      (freePAModule,PathAlgebra,ZZ,HashTable),
	      (freePAModule,PathAlgebraQuotient,ZZ,HashTable),
	      isUniform,
	      (isUniform,PAElement),
	      getSyzygies,
	      (getSyzygies,List,List),
	      interreduce, 
	      (interreduce,List),
	      removeZeroes,
	      reduceOverlap,
	      (reduceOverlap,Sequence,List,List,List),
	      getSyzygiesMatrix,
	      getMinSyzygies,
	      (getMinSyzygies,List,List),
	      getMinSyzMatrix,
              (getMinSyzMatrix,PAMatrix,List),
	      (isSubword,PAElement,PAModMon),
	      (paMatrix,List,HashTable)
	      }	     

beginDocumentation()

--------------------
------ Types -------
--------------------

doc ///
  Key
    PathAlgebras
  Headline
    Data types and basic functions on path algebras.
  Description
    Text
      This package is used to define and manipulate Path Algebras.
  Subnodes
    "Basic operations on Path Algebras"

///

doc ///
  Key
    "Basic operations on Path Algebras"
  Headline
    Outlines some basic operations on Path Algebras
  Description
    Text
     "TO DO"
    Example
     "TO DO"
///



doc ///
  Key
    PathAlgebra
  Headline
    The class of all PathAlgebras
  Description
    Text
     Some common ways to create PathAlgebra include @ TO putInPathAlgebra @
  SeeAlso
     "to do"
///

doc ///
  Key
    PAPath
  Headline
    Type of path in a Path Algebra
  Description
    Text
     "to do"
  SeeAlso
     "to do"
///

doc ///
  Key
    PAGraph
  Headline
    Type of a map to or from a Path Algebra
  Description
    Text
     "to do"
  SeeAlso
     "to do"
///

doc ///
  Key
    PAMatrix
  Headline
    Type of a matrix of Path Algebra elements
  Description
    Text
       Common operations on matrices:
    Code
       UL {TO (symbol +, PAMatrix,PAMatrix),
           TO (symbol -, PAMatrix,PAMatrix),
           TO (symbol *, PAMatrix,PAMatrix),
           TO (symbol *, PAMatrix,QQ),
	   TO (symbol *, PAMatrix,ZZ),
	   TO (symbol *, QQ,PAMatrix),
	   TO (symbol *, ZZ,PAMatrix),
	   TO (symbol ^, PAMatrix,ZZ),
	   }
///


doc ///
  Key
    PAMap
  Headline
    Type of a map to and from a Path Algebra.
  Description
    Text
       A map F:R->S where R and S are Path Algebras. The target map is given first. 
    Text
       Common ways to make (and use) an PAMap include
    Code
       UL {TO (paMap, PathAlgebra,PathAlgebra,List,List)}
///

	
doc ///
  Key
    PAIdeal
  Headline
    Type of a Path Algebra ideal
  Description
    Text
     "to do"
  SeeAlso
     "to do"
///

doc ///
  Key
    PAElement
  Headline
    Type of an element in a Path Algebra
  Description
    Text
     "..."
  SeeAlso
     "to do"
///

doc ///
  Key
    PAModMon
  Headline
    Type of an PAModule
  Description
    Text
     "..."
  SeeAlso
     "to do"
///


doc ///
  Key
    PathAlgebraQuotient
  Headline
    Type of a Path Algebra
  Description
    Text 
     "..."  
  SeeAlso
     "to do" 
///

doc ///
  Key
    PAVector
  Headline
    Type of a Path Algebra Vector
  Description
    Text 
       Common operations on vectors:
    Code
       UL {TO (symbol ==, PAVector,PAVector),
           TO (symbol ?, PAVector,PAVector)
	   }
///

doc ///
  Key
    PAModule
  Headline
    Type of a Path Algebra Module
  Description
    Text
     "..."
  SeeAlso
     "to do"
///

doc ///
   Key
     (symbol +, PAMatrix, PAMatrix)
   Headline
     Add PAMatrices
   Usage
     L = M + N
   Inputs
     M : PAMatrix
     N : PAMatrix
   Outputs
     L : PAMatrix
   Description
    Text
       This adds PAMatrices.
    Example
       M = matrix {{1,0},{1,1}}
       G = paGraph({v,w},{e,f,g},M)
       R = QQ
       A = R G
       N1 = paMatrix{{e,g}}
       N2 = paMatrix{{f,g}}
       N1 + N2
///

doc ///
   Key
     (symbol -, PAMatrix, PAMatrix)
   Headline
     Subtract PAMatrices
   Usage
     L = M - N
   Inputs
     M : PAMatrix
     N : PAMatrix
   Outputs
     L : PAMatrix
   Description
    Text
       This subtracts PAMatrices.
    Example
       M = matrix {{1,0},{1,1}}
       G = paGraph({v,w},{e,f,g},M)
       R = QQ
       A = R G
       N1 = paMatrix{{e,g}}
       N2 = paMatrix{{f,g}}
       N1 - N2
///

doc ///
   Key
     (symbol *, PAMatrix, PAMatrix)
   Headline
     Product of PAMatrices
   Usage
     L = M * N
   Inputs
     M : PAMatrix
     N : PAMatrix
   Outputs
     L : PAMatrix
   Description
    Text
       This command allows for the product of composable PAMatrices (or ordinary matrices over the base).
    Example
       M = matrix {{1,0},{1,1}}
       G = paGraph({v,w},{e,f,g},M)
       R = QQ
       A = R G
       N1 = paMatrix{{e,g}}
       N2 = paMatrix{{f,g}}
       N1 * (transpose N2)
///

doc ///
   Key
     (symbol *, PAMatrix, QQ)
   Headline
      Product of PAMatrices
   Usage
     L = M * c
   Inputs
     M : PAMatrix
     c : QQ
   Outputs
     L : PAMatrix
   Description
    Text
       This command allows for the scaling of an @ TO PAMatrix @ by an element in @ TO QQ @.
    Example
       M = matrix {{1,0},{1,1}}
       G = paGraph({v,w},{e,f,g},M)
       R = QQ
       A = R G
       M = paMatrix {{e^2,e*f},{0_A,-g^2}} 
       M*(2/3)
///

doc ///
   Key
     (symbol *, QQ, PAMatrix)
   Headline
      Product of PAMatrices
   Usage
     L = c * M
   Inputs
     c : QQ
     M : PAMatrix
   Outputs
     L : PAMatrix
   Description
    Text
       This command allows for the scaling of an @ TO PAMatrix @ by an element in @ TO QQ @.
    Example
       M = matrix {{1,0},{1,1}}
       G = paGraph({v,w},{e,f,g},M)
       R = QQ
       A = R G
       M = paMatrix {{e^2,e*f},{0_A,-g^2}} 
       (2/3)*M
///

doc ///
   Key
     (symbol *, ZZ, PAMatrix)
   Headline
      Product of PAMatrices
   Usage
     L = c*M
   Inputs
     c : ZZ
     M : PAMatrix
   Outputs
     L : PAMatrix
   Description
    Text
       This command allows for the scaling of an @ TO PAMatrix @ by an element in @ TO ZZ @.
    Example
       M = matrix {{1,0},{1,1}}
       G = paGraph({v,w},{e,f,g},M)
       R = QQ
       A = R G
       M = paMatrix {{e^2,e*f},{0_A,-g^2}} 
       2*M
///
doc ///
   Key
     (symbol *, PAMatrix,ZZ)
   Headline
      Product of PAMatrices
   Usage
     L = M*c
   Inputs
     M : PAMatrix
     c : ZZ
   Outputs
     L : PAMatrix
   Description
    Text
       This command allows for the scaling of an @ TO PAMatrix @ by an element in @ TO ZZ @.
    Example
       M = matrix {{1,0},{1,1}}
       G = paGraph({v,w},{e,f,g},M)
       R = QQ
       A = R G
       M = paMatrix {{e^2,e*f},{0_A,-g^2}} 
       M*2
///

doc ///
   Key
      (symbol ^, PAMatrix, ZZ)
   Headline
      Exponentiate an PAMatrix
   Usage
      L = M^n
   Inputs
      M : PAMatrix
      n : ZZ
   Outputs
      L : PAMatrix
   Description
      Text
         This exponentiates an NCMatrix.
	 The input is assumed to be a nonnegative integer at this time.
      Example
         M = matrix {{1,0},{1,1}}
	 G = paGraph({v,w},{e,f,g},M)
	 R = QQ
	 A = R G
	 M2 = paMatrix {{1_A,f},{0_A,g}}
	 M2
	 M2^2
	 M2^3
///

doc ///
   Key
      (transpose, PAMatrix)
   Headline
      Transposes an PAMatrix
   Usage
      L = transpose M
   Inputs
      M : PAMatrix
   Outputs
      L : PAMatrix
   Description
      Text
         This command transposes an PAMatrix
      Example
         M = matrix {{1,0},{1,1}}
	 G = paGraph({v,w},{e,f,g},M)
	 R = QQ
	 A = R G
	 M2 = paMatrix {{1_A,f},{0_A,g}}
	 M2
	 transpose M2
	 M2^2
	 transpose M2^2
	 M2^3
	 transpose M2^3
///

doc ///
   Key
      (symbol -, PAMatrix)
   Headline
      Negates PAMatrices
   Usage
     L = -M
   Inputs
     M : PAMatrix
   Outputs
     L : PAMatrix
   Description
      Text
         This negates PAMatrices.
      Example
         M = matrix {{1,0},{1,1}}
       	 G = paGraph({v,w},{e,f,g},M)
       	 R = QQ
       	 A = R G
       	 N1 = paMatrix{{e,g}}
    	 -N1
///

doc ///
   Key
      (entries, PAMatrix)
   Headline
      Returns the entries of the PAMatrix
   Usage
      L = entries M
   Inputs
      M : PAMatrix
   Outputs
      L : List
   Description
      Text
         Returns the entries of the PAMatrix as a doubly nested list.
      Example
         M = matrix {{1,0},{1,1}}
	 G = paGraph({v,w},{e,f,g},M)
	 R = QQ
	 A = R G
	 N = paMatrix {{1_A,f},{0_A,g}}
	 entries N
///

doc ///
   Key
      (symbol _, PAMatrix, Sequence)
   Headline
      Select some columns of an PAMatrix
   Usage
      L = M_cols
   Inputs
      M : PAMatrix
      cols : Sequence
   Outputs
      L : PAMatrix
   Description
      Text
         This command selects some columns of an PAMatrix.
      Example
         M = matrix {{1,0},{1,1}}
	 G = paGraph({v,w},{e,f,g},M)
	 R = QQ
	 A = R G
	 N = paMatrix {{1_A,f},{0_A,g}}
	 N_(1,1)
///

doc ///
   Key
      (symbol |, PAMatrix, PAMatrix)
   Headline
      join PAMatrices horizontally
   Usage
      L = M|N
   Inputs
      M : PAMatrix
      N : PAMatrix
   Outputs
      L : PAMatrix
   Description
      Text
         This command join PAMatrices horizontally
      Example
         M = matrix {{1,0},{1,1}}
	 G = paGraph({v,w},{e,f,g},M)
	 R = QQ
	 A = R G
	 M = paMatrix {{e,f},{e,g}}
	 N = paMatrix {{1_A,f},{0_A,g}}
	 L = M|N
///

doc ///
   Key
      (symbol ||, PAMatrix, PAMatrix)
   Headline
      join PAMatrices horizontally
   Usage
      L = M||N
   Inputs
      M : PAMatrix
      N : PAMatrix
   Outputs
      L : PAMatrix
   Description
      Text
         This command join PAMatrices virtically
      Example
         M = matrix {{1,0},{1,1}}
	 G = paGraph({v,w},{e,f,g},M)
	 R = QQ
	 A = R G
	 M = paMatrix {{e,f},{e,g}}
	 N = paMatrix {{1_A,f},{0_A,g}}
	 L = M||N
///


doc ///
   Key
     (symbol +, PAIdeal, PAIdeal)
   Headline
     Add PAIdeal
   Usage
     L = M + N
   Inputs
     M : PAIdeal
     N : PAIdeal
   Outputs
     L : PAIdeal
   Description
    Text
       This adds PAIdeal.

///



doc ///
   Key
      paMatrix
      (paMatrix,List)
   Headline
      Create an PAMatrix
   Usage
      M = paMatrix entriesList
   Inputs
      entriesList : List
   Outputs
      M : PAMatrix
   Description
      Text
         This command creates an PAMatrix.
      Example
         adj = matrix {{1,0},{1,1}}
         G = paGraph({v,w},{e,f,g},adj)
	 R = QQ
	 A = R G
	 M = paMatrix {{e^2,e*f},{0_A,-g^2}} 
///

doc ///
   Key
      paModMon
      (paModMon,PAModule,ZZ,PAPath)
   Headline
      Create an PA Module Monomial
   Usage
      m = paModMon(M,x,L)
   Inputs
      M:PAModule
      x:ZZ
      L:PAPath
   Outputs
      m : PAModMon
   Description
      Text
         This command creates an PAModMon.
      Example
         adj = matrix {{1,0},{1,1}}
         G = paGraph({v,w},{e,f,g},adj)
	 R = QQ
	 A = R G
         L = leadPath(e*f*g)
	 x = 0
	 M = A^1
	 m = paModMon(M,x,L)
///

--------------------
---- Functions -----
--------------------

doc ///
  Key
    paBasis
    (paBasis,ZZ,PathAlgebra)
    [paBasis,Strategy]
  Headline
    Get a basis for a particular degree or length of a Path Algebra.
  Usage
    L = paBasis(n,A,Option)
  Inputs
    n:ZZ
    A:PathAlgebra
  Outputs
    L:List
      The basis of the desired degree or length of the Path Algebra.
  Description
    Text 
      There are two Strategy that are available for this function, "Degree" and "Length".
      The "Degree" Option is to allow for the retrieval of a basis of a particular degree of a @ TO PathAlgebra @.
      The "Length" Option is to allow for the retrieval of a basis of a particular length of a @ TO PathAlgebra @.
           
    Example      
      M = matrix {{1,0},{1,1}}
      G = paGraph({v,w},{e,f,g},M,Weights=>{2,1,1},Degrees=>{2,1,1})
      R = QQ
      A = R G
      paBasis(1,A,Strategy=>"Degree")
      paBasis(2,A,Strategy=>"Degree")
      paBasis(7,A,Strategy=>"Degree")
      paBasis(1,A,Strategy=>"Length")
      paBasis(7,A,Strategy=>"Length")

///

doc ///
  Key
    (paBasis,ZZ,PathAlgebra,List)
  Headline
    Get a basis for a particular degree or length of a quotient of the Path Algebra.
  Usage
    L = paBasis(n,A,M,Option)
  Inputs
    n:ZZ
    A:PathAlgebra
    M:List
      The list of PAElements whose lead terms the basis avoids.
  Outputs
    L:List
      The basis of the desired degree or length of the Path Algebra.
  Description
    Text
      This gives a basis of the quotient R/M. Usually, this is used when M is the List of lead terms of some other ideal I.
    Text 
      There are two Strategy that are available for this function, "Degree" and "Length".
      The "Degree" Option is to allow for the retrieval of a basis of a particular degree of a @ TO PathAlgebra @.
      The "Length" Option is to allow for the retrieval of a basis of a particular length of a @ TO PathAlgebra @.
           
    Example      
      N = matrix {{1,0},{1,1}}
      G = paGraph({v,w},{e,f,g},N,Weights=>{2,1,1},Degrees=>{2,1,1})
      R = QQ
      A = R G
      paBasis(7,A,{e},Strategy=>"Length")
      paBasis(7,A,{e*f},Strategy=>"Length")
      paBasis(7,A,{e*f},Strategy=>"Degree")
      paBasis(8,A,{e},Strategy=>"Degree")
      paBasis(8,A,{e*f},Strategy=>"Degree")

///


doc ///
  Key
    paMap
    (paMap, PathAlgebra,PathAlgebra,List,List)
  Headline
    Get Path Algebra mapping from one PathAlgebra to another
  Usage
    L = paMap(A,B,M,N)
  Inputs
    A:PathAlgebra
    B:PathAlgebra
    M:List
    N:List
  Outputs
    L:PAMap

///

doc ///
  Key
    (paMap,PathAlgebra,PathAlgebraQuotient,List,List)
  Headline
    Get Path Algebra mapping from a PathAlgebraQuotient to a PathAlgebra
  Usage
    L = paMap(A,B,M,N)
  Inputs
    A:PathAlgebra
    B:PathAlgebraQuotient
    M:List
    N:List
  Outputs
    L:PAMap

///


doc ///
  Key
    (paMap,PathAlgebraQuotient,PathAlgebra,List,List)
  Headline
    Get Path Algebra mapping from a PathAlgebra to a PathAlgebraQuotient
  Usage
    L = paMap(A,B,M,N)
  Inputs
    A:PathAlgebraQuotient
    B:PathAlgebra
    M:List
    N:List
  Outputs
    L:PAMap

///


doc ///
  Key
    (paMap,PathAlgebraQuotient,PathAlgebraQuotient,List,List)
  Headline
    Get Path Algebra mapping from one PathAlgebraQuotient to another
  Usage
    L = paMap(A,B,M,N)
  Inputs
    A:PathAlgebraQuotient
    B:PathAlgebraQuotient
    M:List
    N:List
  Outputs
    L:PAMap

///



doc ///
  Key
    isSubword
    (isSubword,List,List)
  Headline
    Determines if a list of edges is a subword of another.
  Usage
    isSub = isSubword(p,q)
  Inputs
    p:List
    q:List
  Outputs
    isSub:Sequence
          (Boolean,ZZ)
  Description
    Text
      This function determines if the first list of edges is a subword of the second. p and q are the edgelist of PAElements.
    Text
      Return value is (Boolean,ZZ) where the boolean is true if p is a subword of q and false otherwise,
      and the integer is the first occurrence (from the left) of p in q. Returns -1 if not subword.
    Example
      p={0,0,0,1,2}
      q={0,1,2}
      isSubword(p,q)
    Example
      p={0,1,2}
      q={0,0,0,1,2}
      isSubword(p,q)
///


doc ///
  Key
    paBasisLength
    (paBasisLength,ZZ,PathAlgebra)
  Headline
    Get a basis for a particular length of a Path Algebra.
  Usage
    L = paBasisLength(n,A)
  Inputs
    n:ZZ
    A:PathAlgebra
  Outputs
    L:List
      The basis of the desired length of the Path Algebra.
  Description
    Text 
      This function is to allow for the retrieval of a basis of a particular length of a @ TO PathAlgebra @.
           
    Example      
      M = matrix {{1,0},{1,1}}
      G = paGraph({v,w},{e,f,g},M,Weights=>{2,1,1},Degrees=>{2,1,1})
      R = QQ
      A = R G
      paBasisLength(1,A)
      paBasisLength(7,A)


///

doc ///
  Key
    (paBasisLength,ZZ,PathAlgebra,List)
  Headline
    Get a basis for a particular length of a quotient of Path Algebra.
  Usage
    L = paBasisLength(n,A,M)
  Inputs
    n:ZZ
    A:PathAlgebra
    M:List
      A List of PAElements whose lead terms the basis avoids.
  Outputs
    L:List
      The basis of the desired length of the Path Algebra.
  Description
    Text
      This function is to allow for the retrieval of a basis of a particular length of the quotient R/M of a @ TO PathAlgebra @.
      Usually, this is used when M is the list of lead terms of some other ideal I.
           
    Example      
      N = matrix {{1,0},{1,1}}
      G = paGraph({v,w},{e,f,g},N,Weights=>{2,1,1},Degrees=>{2,1,1})
      R = QQ
      A = R G
      paBasisLength(1,A,{e})
      paBasisLength(7,A,{e})
      paBasisLength(7,A,{e*f})

///


doc ///
  Key
    paBasisDegree
    (paBasisDegree,ZZ,PathAlgebra)
  Headline
    Get a basis for a particular degree of a Path Algebra.
  Usage
    L = paBasisDegree(n,A)
  Inputs
    n:ZZ
    A:PathAlgebra
  Outputs
    L:List
      The basis of the desired degree of the Path Algebra.
  Description
    Text 
      This function is to allow for the retrieval of a basis of a particular degree of a @ TO PathAlgebra @.
           
    Example      
      N = matrix {{1,0},{1,1}}
      G = paGraph({v,w},{e,f,g},N,Weights=>{2,1,1},Degrees=>{2,1,1})
      R = QQ
      A = R G
      paBasisDegree(1,A)
      paBasisDegree(2,A)
      paBasisDegree(7,A)


///



doc ///
  Key
    (paBasisDegree,ZZ,PathAlgebra,List)
  Headline
    Get a basis for a particular degree of a quotient of Path Algebra.
  Usage
    L = paBasisDegree(n,A,M)
  Inputs
    n:ZZ
    A:PathAlgebra
    M:List
      A List of PAElements whose lead terms the basis avoids.
  Outputs
    L:List
      The basis of the desired degree of the Path Algebra.
  Description
    Text
      This function is to allow for the retrieval of a basis of a particular degree of the quotient R/M of a @ TO PathAlgebra @.
      Usually, this is used when M is the list of lead terms of some other ideal I.
           
    Example      
      N = matrix {{1,0},{1,1}}
      G = paGraph({v,w},{e,f,g},N,Weights=>{2,1,1},Degrees=>{2,1,1})
      R = QQ
      A = R G
      paBasisDegree(1,A,{e})
      paBasisDegree(7,A,{e})
      paBasisDegree(7,A,{e*f})

///


doc ///
  Key
    isSubpath
    (isSubpath,PAElement,PAElement)
  Headline
    Determines if a PAElement is a subpath of another.
  Usage
    isSub = isSubpath(p,q)
  Inputs
    p:PAElement
    q:PAElement
  Outputs
    isSub:Sequence

  Description
    Text
      This function determines if a PAElement is a subpath of another. p and q are PAElements.
    Text
      Return value is (Boolean,symbol,PAElement,PAElement) where the boolean is true if p is a subpath of q or q is a subpath of p, and false otherwise.
      And the symbol determines which PAElement is a subpath: returns > if p is a subpath of q and returns < if q is a subpath of p, and returns "not Subpath" if not both.
      The returned PAElements are u and v such that p=u*q*v or q=u*p*v.
    Example
      M = matrix {{1,0},{1,1}}
      G = paGraph({v,w},{e,f,g},M)
      R = QQ
      A = R G
      p=e*e*f*g
      q=e*f
      isSubpath(p,q)
      isSubpath(q,p)
      m=e*e*e*f*g
      n=e*f*g*g*g
      isSubpath(m,n)

///


doc ///
  Key
    isSubpathOnly
    (isSubpathOnly,PAElement,PAElement)
  Headline
    Determines if a PAElement is a subpath of another.
  Usage
    isSub = isSubpathOnly(p,q)
  Inputs
    p:PAElement
    q:PAElement
  Outputs
    isSub:Boolean
  Description
    Text
      This function determines if a PAElement is a subpath of another. p and q are PAElements.
    Text
      Return value is Boolean where it is true if p is a subpath of q or q is a subpath of p. Return false if not subpath.

    Example
      M = matrix {{1,0},{1,1}}
      G = paGraph({v,w},{e,f,g},M)
      R = QQ
      A = R G
      p=e*e*f*g
      q=e*f
      isSubpathOnly(p,q)
      isSubpathOnly(q,p)
      m=e*e*e*f*g
      n=e*f*g*g*g
      isSubpathOnly(m,n)
///

doc ///
  Key
    weight
    (weight,PAElement)
  Headline
    Calculate the weight of a PAElement.
  Usage
    w = weight(p)
  Inputs
    p:PAElement
  Outputs
    w:ZZ
  Description
    Text
      This function calculates the weight of a PAElement. Return value is the sum of weight of each edges of this PAElement.

    Example
      M = matrix {{1,0},{1,1}}
      G = paGraph({v,w},{e,f,g},M,Weights=>{2,1,1},Degrees=>{2,1,1})
      R = QQ
      A = R G
      p=e*e*f*g
      weight(p)

///

doc ///
  Key
    (weight,PAVector)
  Headline
    Calculate the weight of a PAVector.
  Usage
    w = weight(p)
  Inputs
    p:PAVector
  Outputs
    w:ZZ
  Description
    Text
      This function calculates the weight of a PAVector. Return value is the sum of weight of each edges of this PAVector.

///

doc ///
  Key
   leadTerm
    (leadTerm,PAElement)
  Headline
    Find the lead term of a PAElement
  Usage
    L = leadTerm(f)
  Inputs
    f:PAElement
  Outputs
    L:PAElement
  Description
    Text
      This function finds the lead term of a PAElment.

    Example
      M = matrix {{1,0},{1,1}}
      G = paGraph({v,w},{e,f,g},M,Weights=>{2,1,1},Degrees=>{2,1,1})
      R = QQ
      A = R G
      f=e*e*f*g+e*f*g+f*g
      L=leadTerm(f)

///

doc ///
  Key
    (leadTerm,PAVector)
  Headline
    Find the lead term of a PAVector
  Usage
    L = leadTerm(f)
  Inputs
    f:PAVector
  Outputs
    L:PAVector
  Description
    Text
      This function finds the lead term of a PAVector.

    Example
      "TO DO"

///


doc ///
  Key
   leadPath
    (leadPath,PAElement)
  Headline
    Find the lead path of a PAElement
  Usage
    L = leadPath(f)
  Inputs
    f:PAElement
  Outputs
    L:PAPath
  Description
    Text
      This function finds the lead path of a PAElment.
    Example
      "TO DO"


///


doc ///
  Key
   leadEdgeList
    (leadEdgeList,PAElement)
  Headline
    Find the lead edgelist of a PAElement
  Usage
    L = leadEdgeList(f)
  Inputs
    f:PAElement
  Outputs
    L:List
  Description
    Text
      This function finds the lead edgelist of a PAElement.
    Example
      "Need Examples here"


///

doc ///
  Key
   leadTermCoeff
    (leadTermCoeff,PAElement)
  Headline
    Find the coefficient of the lead term of a PAElement.
  Usage
    L = leadTermCoeff(f)
  Inputs
    f:PAElement
  Outputs
    L:RingElement
  Description
    Text
      This function finds the coefficient of the lead term of a PAElement.

    Example
      M = matrix {{1,0},{1,1}}
      G = paGraph({v,w},{e,f,g},M,Weights=>{2,1,1},Degrees=>{2,1,1})
      R = QQ
      A = R G
      f=2*e*e*f*g+e*f*g+f*g
      leadTermCoeff(f)
///


doc ///
  Key
    (leadTermCoeff,PAVector)
  Headline
    Find the coefficient of the lead term of a PAVector.
  Usage
    L = leadTermCoeff(f)
  Inputs
    f:PAVector
  Outputs
    L:RingElement
  Description
    Text
      This function finds the coefficient of the lead term of a PAVector.

    Example
      M = matrix {{1,0},{1,1}}
      G = paGraph({v,w},{e,f,g},M,Weights=>{2,1,1},Degrees=>{2,1,1})
      R = QQ
      A = R G
      B = A^5
      e0 = B_0
      L = 5*e0*e*f*g+4*e0*e*f+3*e0*f*g
      leadTermCoeff(L)
///	      
	      
	      
doc ///
  Key
    pathDegree
    (pathDegree,PAElement)
  Headline
    Find the degree of the lead term of a PAElement.
  Usage
    L = pathDegree(f)
  Inputs
    f:PAElement
  Outputs
    L:ZZ
  Description
    Text
      This function finds the degree of the lead term of a PAElement.
    Example
      "TO DO"
///

doc ///
  Key
    (pathDegree,PAVector)
  Headline
    Find the degree of the lead term of a PAVector.
  Usage
    L = pathDegree(f)
  Inputs
    f:PAVector
  Outputs
    L:ZZ
  Description
    Text
      This function finds the degree of the lead term of a PAVector.
    Example
      "TO DO"
///


doc ///
  Key
    paGraph
    (paGraph,List,List,Matrix)
    [paGraph,Degrees]
    [paGraph,Weights]
  Headline
   Define a Path Algebra Graph.
  Usage
    L = paGraph(verts,edges,adj)
  Inputs
    verts:List
    edges:List
    adj:Matrix
  Outputs
    L:PAGraph
  Description
    Text 
      "to do"        
    Example      
      M = matrix {{1,0},{1,1}}
      G = paGraph({v,w},{e,f,g},M,Weights=>{2,1,1},Degrees=>{2,1,1})


///


doc ///
  Key
    paIdeal
    (paIdeal,List)
  Headline
   Generates a Ideal of Path Algebras
  Usage
    L = paIdeal(l)
  Inputs
    l:List
  Outputs
    L:PAIdeal
  Description
    Text 
      "to do"        
    Example      
      M = matrix {{3}}
      G = paGraph({v},{a,b,c},M,Weights=>{1,1,1})
      kk = ZZ/32003
      R = kk G
      I = paIdeal {2*a*b + 3*b*a + 5*c^2,
             2*b*c + 3*c*b + 5*a^2,
             2*c*a + 3*a*c + 5*b^2}
	 
///

doc ///
  Key
    (length,PAPath)
  Headline
    Outputs the length of a PAPath.
  Usage
    L=length(f)
  Inputs
    f:PAPath
  Outputs
    L:ZZ
///

doc ///
  Key
    (length,PAModMon)
  Headline
    Outputs the length of a PAModMon.
  Usage
    L=length(f)
  Inputs
    f:PAModMon
  Outputs
    L:ZZ
///

doc ///
  Key
    (length,PAElement)
  Headline
    Outputs the length of a PAElement.
  Usage
    L=length(f)
  Inputs
    f:PAElement
  Outputs
    L:ZZ
///


doc ///
  Key
    (net,PAGraph)
  Headline
    Outputs the pertinent information about a PAGraph
  Usage
    net B
  Inputs
    B:PAGraph
///

doc ///
  Key
    (net,PAPath)
  Headline
    Outputs the pertinent information about a PAPath
  Usage
    net C
  Inputs
    C:PAPath
///

doc ///
  Key
    (net,PAElement)
  Headline
    Outputs the pertinent information about a PAElement
  Usage
    net D
  Inputs
    D:PAElement
///

doc ///
  Key
    (net,PAMatrix)
  Headline
    Outputs the pertinent information about a PAMatrix
  Usage
    net E
  Inputs
    E:PAMatrix
///

doc ///
  Key
    (net,PAVector)
  Headline
    Outputs the pertinent information about a PAVector
  Usage
    net F
  Inputs
    F:PAVector
///

doc ///
  Key
    (net,PAModule)
  Headline
    Outputs the pertinent information about a PAModule
  Usage
    net G
  Inputs
    G:PAModule
///

doc ///
  Key
    (ring,PAElement)
  Headline
    Outputs the ring information about a PAElement
  Usage
    R=ring f
  Inputs
    f:PAElement
  Outputs
    R:Ring
///

doc ///
  Key
    buchAlgorithm
    (buchAlgorithm,List)
  Headline
    Computes the Grobner Basis
  Usage
    G=buchAlgorithm(L)
  Inputs
    L:List
  Outputs
    G:List
  Description
    Text
     "to do"
    Example
     M = matrix {{1,0},{1,1}}
     G = paGraph({v,w},{e,f,g},M,Weights=>{2,1,1},Degrees=>{2,1,1})
     R = QQ
     A = R G
     buchAlgorithm({e*f*g+f,e*e*f})

///

doc ///
  Key
    (buchAlgorithm,PAIdeal)
    [buchAlgorithm,DegreeLimit]
  Headline
    Computes the Grobner Basis
  Usage
    G=buchAlgorithm(I)
  Inputs
    L:PAIdeal
  Outputs
    G:PAIdeal
  Description
    Text
     "to do"
    Example
     M = matrix {{3}}
     G = paGraph({v},{a,b,c},M,Weights=>{1,1,1})
     kk = ZZ/32003
     R = kk G
     I = paIdeal {2*a*b + 3*b*a + 5*c^2,
             2*b*c + 3*c*b + 5*a^2,
             2*c*a + 3*a*c + 5*b^2}
     gbI = netList sort buchAlgorithm(I,DegreeLimit=>7)
     
///


doc ///
  Key
    (symbol ==, PAElement,PAElement)
    (symbol ==, PAElement,ZZ)
    (symbol ==, ZZ,PAElement)
  Headline
    Basic compares with PAElements.
  Usage
    b=(f==g)
  Inputs
    f : PAElement
    g : PAElement
  Outputs
    b:Boolean
  Description
    Text
      Compares two PAElements. Returns true if two PAElements are the same, false otherwise.
      
///


doc ///
  Key
    (symbol ?, PAElement,PAElement)
  Headline
    Basic compares with PAElements.
  Usage
    b=(f?g)
  Inputs
    f : PAElement
    g : PAElement
  Outputs
    b:Boolean
  Description
    Text
      Compares two PAElements based on term orders.    
///

doc ///
  Key
    (symbol ?, PAVector,PAVector)
  Headline
    Basic compares with PAVector.
  Usage
    b=(f?g)
  Inputs
    f : PAVector
    g : PAVector
  Outputs
    b : Boolean
  Description
    Text
      Compares two PAVectors based on term orders.    
///

doc ///
  Key
    (symbol ==, PAVector,PAVector)
    (symbol ==, PAVector,ZZ)
  Headline
    Basic compares with PAVectors.
  Usage
    b=(f==g)
  Inputs
    f : PAVector
    g : PAVector
  Outputs
    b : Boolean
  Description
    Text
      Compares two PAVectors. Returns true if two PAVectors are the same, false otherwise.
      
///

doc ///
  Key
    (symbol ?, PAModMon,PAModMon)
  Headline
    Basic compares with PAModMons.
  Usage
    b=(f?g)
  Inputs
    f : PAModMon
    g : PAModMon
  Outputs
    b : Boolean
  Description
    Text
      Compares two PAModMons based on term orders.    
///

doc ///
  Key
    (symbol ==, PAModMon,PAModMon)
  Headline
    Basic compares with Path Algebra Module Monomials.
  Usage
    b=(f==g)
  Inputs
    f : PAModMon
    g : PAModMon
  Outputs
    b : Boolean
  Description
    Text
      Compares two PAModMons. Returns true if two PAElements are the same, false otherwise.
      
///
doc ///
  Key
    (symbol ==, PAPath,PAPath)
  Headline
    Basic compares with PAPaths.
  Usage
    b=(f==g)
  Inputs
    f : PAPath
    g : PAPath
  Outputs
    b:Boolean
  Description
    Text
      Compares two PAElements. Returns true if two PAElements are the same, false otherwise.
      
///

doc ///
  Key
    (symbol ?, PAPath,PAPath)
  Headline
    Basic compares with PAPaths.
  Usage
    b=(f?g)
  Inputs
    f : PAPath
    g : PAPath
  Outputs
    b:Boolean
  Description
    Text
      Compares two PAPaths based on term orders.    
///

doc ///
  Key
    (symbol %,PAElement,List)
  Headline
    Find the remainder of a PAElement divided by a List of PAElements.
  Usage
    r=f%L
  Inputs
    f:PAElement
    g:List
        A List of PAElements
  Outputs
    r:PAElement   
///

doc ///
  Key
    (symbol %,PAElement,PAIdeal)
  Headline
    Find the remainder of a PAElement mod a PAIdeal
  Usage
    r=f%L
  Inputs
    f:PAElement
    g:PAIdeal
  Outputs
    r:PAElement   
///


doc ///
   Key
     (symbol /, PathAlgebra, PAIdeal)
   Headline
     Construct a PathAlgebraQuotient
   Usage
     B = A/I
   Inputs
     A : PathAlgebra
     I : PAIdeal
   Outputs
     B : PathAlgebraQuotient
   Description
      Text
         This is one way to create a quotient of the Path Algebra modulo some PAIdeals.
      Example
         M = matrix {{3}}
	 G = paGraph({v},{a,b,c},M,Weights=>{1,1,1})
	 kk = ZZ/32003
	 R = kk G
	 I = paIdeal {2*a*b + 3*b*a + 5*c^2,
             2*b*c + 3*c*b + 5*a^2,
             2*c*a + 3*a*c + 5*b^2}
	 S = R/I
///


doc ///
   Key
     areComposable
     (areComposable,PAGraph, PAPath, PAPath)
   Headline
     Determine whether two PAPath are composable
   Usage
     AC=areComposable(G,p,q)
   Inputs
     G:PAGraph
     p:PAPath
     q:PAPath
   Outputs
     AC: Boolean
   Description
      Text
        This function returns true or false, depending on whether p and q are composable
      Example
        M = matrix {{1,0},{1,1}}
	G = paGraph({v,w},{e,f,g},M,Weights=>{2,1,1},Degrees=>{2,1,1})
	R = QQ
	A = R G
        areComposable(G,leadPath(e*f*g),leadPath(g*g*g)) 
        areComposable(G,leadPath(g*g*g),leadPath(e*f*g))  
///

doc ///
   Key
     composePath
     (composePath,PAPath, PAPath)
   Headline
     Compose two PAPaths
   Usage
     comp=composePath(p,q)
   Inputs
     p:PAPath
     q:PAPath
   Outputs
     comp:PAPath
   Description
      Text 
       This function assumes p and q are composable paths and produces the composition of p and q.
      Example
        M = matrix {{1,0},{1,1}}
	G = paGraph({v,w},{e,f,g},M,Weights=>{2,1,1},Degrees=>{2,1,1})
	R = QQ
	A = R G
	--composePath(leadPath(e*f*g),leadPath(g*g*g))
	L=paPath(G,{0,1,2})
	J=paPath(G,{2,2,2})
	composePath(L,J)
///

doc ///
  Key
    rightLexCompare
    (rightLexCompare,PAPath,PAPath)
  Headline
    Right length-lexicographic order
  Usage
    compare=rightLexCompare(p,q)
  Inputs
    p:PAPath
    q:PAPath
  Outputs
    compare:Boolean
  Description
    Text 
       "To Do"       
      
///


doc ///
  Key
    leftLexCompare
    (leftLexCompare,PAPath,PAPath)
  Headline
    Left length-lexicographic order
  Usage
    compare=leftLexCompare(p,q)
  Inputs
    p:PAPath
    q:PAPath
  Outputs
    compare:Boolean
  Description
    Text 
        "To Do"     
      
///


doc ///
  Key
    leftWeightLexCompare
    (leftWeightLexCompare,List,PAPath,PAPath)
  Headline
    Left weight-lexicographic order
  Usage
    compare=leftWeightLexCompare(w,p,q)
  Inputs
    w:List
    p:PAPath
    q:PAPath
  Outputs
    compare:Boolean
  Description
    Text 
       "To Do"      
      
///

doc ///
  Key
    rightWeightLexCompare
    (rightWeightLexCompare,List,PAPath,PAPath)
  Headline
    Right weight-lexicographic order
  Usage
    compare=rightWeightLexCompare(w,p,q)
  Inputs
    w:List
    p:PAPath
    q:PAPath
  Outputs
    compare:Boolean
  Description
    Text 
       "To Do"      
      
///

doc ///
  Key
    leftWeightReverseLexCompare
    (leftWeightReverseLexCompare,List,PAPath,PAPath)
  Headline
    Left weight-reverse-lexicographic order
  Usage
    compare=leftWeightReverseCompare(w,p,q)
  Inputs
    w:List
    p:PAPath
    q:PAPath
  Outputs
    compare:Boolean
  Description
    Text 
       "To Do"         
      
///

doc ///
  Key
    rightWeightReverseLexCompare
    (rightWeightReverseLexCompare,List,PAPath,PAPath)
  Headline
    Right weight-reverse-lexicographic order
  Usage
    compare=rightWeightReverseLexCompare(w,p,q)
  Inputs
    w:List
    p:PAPath
    q:PAPath
  Outputs
    compare:Boolean
  Description
    Text 
       "To Do"        
      
///



doc ///
  Key
    startVertex
    (startVertex,PAPath)
  Headline
    Find the starting vertex of a PAPath
  Usage
    stv=startVertex(p)
  Inputs
    p:PAPath
  Description
    Text 
       "To Do"        
      
///

doc ///
  Key
    (startVertex,PAElement)
  Headline
    Find the starting vertex of a PAElement
  Usage
    stv=startVertex(p)
  Inputs
    p:PAElement
  Description
    Text 
       "To Do"        
      
///

-*
doc ///
  Key
    (startVertex,PAVector)
  Headline
    Find the starting vertex of a PAVector
  Usage
    stv=startVertex(p)
  Inputs
    p:PAVector
  Description
    Text 
       "To Do"        
      
///
*-

doc ///
  Key
    (startVertex,List)
  Headline
    Find the starting vertex of a list of PAPaths,PAElements,or PAVectors.
  Usage
    stv=startVertex(L)
  Inputs
    L:List
  Description
    Text 
       "To Do"        
      
///

doc ///
  Key
    endVertex
    (endVertex,PAPath)
  Headline
    Find the ending vertex of a PAPath
  Usage
    endv=endVertex(p)
  Inputs
    p:PAPath
  Description
    Text 
       "To Do"        
      
///

doc ///
  Key
    (endVertex,PAElement)
  Headline
    Find the ending vertex of a PAElement
  Usage
    endv=endVertex(p)
  Inputs
    p:PAElement
  Description
    Text 
       "To Do"           
///

doc ///
  Key
    (endVertex,PAVector)
  Headline
    Find the ending vertex of a PAVector
  Usage
    endv=endVertex(p)
  Inputs
    p:PAVector
  Description
    Text 
       "To Do"           
///

doc ///
  Key
    (endVertex,List)
  Headline
    Find the ending vertex of a list of PAPaths,PAElements,or PAVectors.
  Usage
    endv=endVertex(p)
  Inputs
    p:List
  Description
    Text 
       "To Do"           
///


doc ///
   Key
      (terms, PAElement)
   Headline
      Returns the terms of an PAElement
   Usage
     t = terms f
   Inputs
     f : PAElement
   Outputs
     t : List
   Description
      Text
         Returns the list of terms that make up the PAElement. 
	 It is a list of PAElements.
      Example
         adj = matrix {{1,0},{1,1}}
	 G = paGraph({v,w},{e,f,g},adj)
	 R = QQ
	 A = R G
	 L = e*f*g+e*f
	 t = terms L
///

doc ///
   Key
      (terms, PAVector)
   Headline
      Returns the terms of an PAVector
   Usage
     t = terms f
   Inputs
     f : PAVector
   Outputs
     t : List
   Description
      Text
         Returns the list of terms that make up the PAVector. 
	 It is a list of PAElements.
      Example
         "TO DO"
///

doc ///
   Key
     hasDuplicate
     (hasDuplicate,List)
   Headline
     Check if repeated edgelist exists.
   Usage
     dup=hasDuplicate(L)
   Inputs
     L:List
   Outputs
     dup:Boolean
   Description
      Text 
       This function takes a List of edgelists and check if repeated edgelist exists.
      Example
       A={0,1,2}
       B={0,1}
       C={0,0,0,1}
       D={0,1}
       E={B,C}
       F={B,C}
       G={A,B}
       dup=hasDuplicate(E)
       dup2=hasDuplicate(F)
       dup3=hasDuplicate(G)       
///

doc ///
   Key
     findOverlaps
     (findOverlaps,List,List)
   Headline
     Find Overlap Relations
   Usage
     ov=findOverlaps(M,N)
   Inputs
     M:List
     N:List
   Outputs
     ov:List
   Description
      Text 
       This function takes two List of edgelists and return value is a list of positions in p where a (proper nontrivial) suffix of p is a (proper nontrivial) prefix of q.
      Example
       M={0,0,0,1,2,2,2,1,0,1,2} 
       N={0,1,2,2,2,1,0,1,2,0,0}  
       ov=findOverlaps(M,N)
       ov2=findOverlaps(N,M)
      Text
       The following example shows that this function does not allow proper overlaps
      Example
       findOverlaps({1,2,3},{2})
       findOverlaps({1,2,3,1},{3,1})
       findOverlaps({1,2,3,1,3},{3,1,3,4})
///

doc ///
   Key
     (findOverlaps,PAElement,PAModMon)
   Headline
     Find Overlap Relations
   Usage
     ov=findOverlaps(M,N)
   Inputs
     M:PAElement
     N:PAModMon
   Outputs
     ov:PAPath
   Description
      Text 
       This function takes one PAElement and one PAModMon and returns their overlap relation. 
      Example
        adj = matrix {{3}}
	G = paGraph({v},{e,f,g},adj)
	R = QQ
	A = R G
	F0 = A^1
	e0 = F0_0
	f1 = e*f
	f2 = e0*e*f*g
	J = leadModMon f2
	findOverlaps(f1,J)

///


doc ///
   Key
     (findOverlaps,PAModMon,PAElement)
   Headline
     Find Overlap Relations
   Usage
     ov=findOverlaps(M,N)
   Inputs
     M:PAModMon
     N:PAElement
   Outputs
     ov:PAPath
   Description
      Text 
        This function takes one PAElement and one PAModMon and returns their overlap relation. 
      Example
        adj = matrix {{3}}
	G = paGraph({v},{e,f,g},adj)
	R = QQ
	A = R G
	F0 = A^1
	e0 = F0_0
	f1 = e0*e*f
	f2 = e*f*g
        L = leadModMon f1
	findOverlaps(L,f2)

///

doc ///
   Key
     (findOverlaps,PAModMon,PAModMon)
   Headline
     Find Overlap Relations
   Usage
     ov=findOverlaps(M,N)
   Inputs
     M:PAModMon
     N:PAModMon
   Outputs
     ov:PAPath
   Description
      Text 
       This function takes two PAModMons and returns their overlap relation. 
      Example
        adj = matrix {{3}}
	G = paGraph({v},{e,f,g},adj)
	R = QQ
	A = R G
	F0 = A^1
	e0 = F0_0
	f1 = e0*e*f
	f2 = e0*e*f*g
        L = leadModMon f1
	J = leadModMon f2
	findOverlaps(L,J)

///


doc ///
   Key
     (findOverlaps,PAModMon,PAPath)
   Headline
     Find Overlap Relations
   Usage
     ov=findOverlaps(M,N)
   Inputs
     M:PAModMon
     N:PAPath
   Outputs
     ov:PAPath
   Description
      Text 
       This function takes a PAModMon and a PAPath and returns their overlap relation. 
      Example
        adj = matrix {{3}}
	G = paGraph({v},{e,f,g},adj)
	R = QQ
	A = R G
	F0 = A^1
	e0 = F0_0
	f1 = e0*e*f
	f2 = e*f*g
        L = leadModMon f1
	J = leadPath f2
	findOverlaps(L,J)

///

doc ///
   Key
     (findOverlaps,PAPath,PAModMon)
   Headline
     Find Overlap Relations
   Usage
     ov=findOverlaps(M,N)
   Inputs
     M:PAPath
     N:PAModMon

   Outputs
     ov:PAPath
   Description
      Text 
       This function takes a PAModMon and a PAPath and returns their overlap relation. 
      Example
        adj = matrix {{3}}
	G = paGraph({v},{e,f,g},adj)
	R = QQ
	A = R G
	F0 = A^1
	e0 = F0_0
	f1 = e*f
	f2 = e0*e*f*g
        L = leadPath f1
	J = leadModMon f2
	findOverlaps(L,J)

///

doc ///
   Key
     leftOverlaps
     (leftOverlaps,PAElement,PAElement)
   Headline
     Find Overlap Relations of two PAElements
   Usage
     ovl=leftOverlaps(f,g)
   Inputs
     f:PAElement
     g:PAElement

   Outputs
     ovl:Sequence
   Description
      Text 
        This function takes two PAElements and returns their overlap relation with multiplication on the left of the second PAElement. 
      Example 
        M = matrix {{3}}
	G = paGraph({v},{a,b,c},M,Weights=>{1,1,1})
	kk = ZZ/32003
	R = kk G
	I = paIdeal {2*a*b + 3*b*a + 5*c^2,2*b*c + 3*c*b + 5*a^2,2*c*a + 3*a*c + 5*b^2}
	ov = leftOverlaps((I.generators)#1,(I.generators)#2)
	peek ov
///


doc ///
   Key
     rightOverlaps
     (rightOverlaps,PAElement,PAElement)
   Headline
     Find Overlap Relations of two PAElements
   Usage
     ovr=rightOverlaps(f,g)
   Inputs
     f:PAElement
     g:PAElement

   Outputs
     ovr:Sequence
   Description
      Text 
       This function takes two PAElements and returns their overlap relation with multiplication on the right of the second PAElement. 
      Example 
        M = matrix {{3}}
	G = paGraph({v},{a,b,c},M,Weights=>{1,1,1})
	kk = ZZ/32003
	R = kk G
	I = paIdeal {2*a*b + 3*b*a + 5*c^2,2*b*c + 3*c*b + 5*a^2,2*c*a + 3*a*c + 5*b^2}
	ov = rightOverlaps((I.generators)#2,(I.generators)#1)
	peek ov
///


doc ///
   Key
     overlaps
     (overlaps,PAElement,PAElement)
   Headline
     Find Overlap Relations of two PAElements
   Usage
     ov=overlaps(f,g)
   Inputs
     f:PAElement
     g:PAElement

   Outputs
     ov:Sequence
   Description
      Text 
       This function takes two PAElements and find all of their overlap relations. 
      Example
        M = matrix {{3}}
	G = paGraph({v},{a,b,c},M,Weights=>{1,1,1})
	kk = ZZ/32003
	R = kk G
	I = paIdeal {2*a*b + 3*b*a + 5*c^2,2*b*c + 3*c*b + 5*a^2,2*c*a + 3*a*c + 5*b^2}
	ov = overlaps((I.generators)#1,(I.generators)#2)
	peek ov

///

doc ///
   Key
     (overlaps,PAVector,PAElement,ZZ)
   Headline
     Find Overlap Relations of a PAVector with a PAElement
   Usage
     ov=overlaps(f,g,n)
   Inputs
     f:PAVector
     g:PAElement
     n:ZZ
   Outputs
     ov:Sequence
   Description
      Text 
       This function takes two PAElements and find all of their overlap relations. 
      Example
        "TODO"

///
doc ///
   Key     
     (overlaps,Sequence,PAElement,ZZ)
   Headline
     Find Overlap Relations of a PAVector with a PAElement, with tracking info
   Usage
     ov=overlaps(f,g,n)
   Inputs
     f:Sequence
     g:PAElement
     n:ZZ
   Outputs
     ov:Sequence
   Description
      Text 
       This function takes two PAElements and find all of their overlap relations. 
      Example
        "TODO"

///


doc ///
   Key
     leadModMon
     (leadModMon,PAVector)
   Headline
     Find the lead Module Monomial of a PAVector
   Usage
     lmm=leadModMon(f)
   Inputs
     f:PAVector
   Outputs
     lmm:PAModMon
   Description
      Text 
       This function gives the lead module monomial of a PAVector 
      Example
        adj = matrix {{1,0},{1,1}}
	G = paGraph({v,w},{e,f,g},adj)
	R = QQ
	A = R G
	F = A^3
	e0 = F_0
	e1 = F_1
	a = e0*e*f
	b = e1*e*f
	c = e0*e*f*g	  
        d=a+b+c
	lmm=leadModMon(d)

///

doc ///
   Key
     allComponent
     (allComponent,PAVector)
   Headline
     Find all of the components of a PAVector
   Usage
     ac=allComponent(f)
   Inputs
     f:PAVector
   Outputs
     ac:List
   Description
      Text 
       This function returns a list of all the components of a PAVector 
      Example
        adj = matrix {{1,0},{1,1}}
	G = paGraph({v,w},{e,f,g},adj)
	R = QQ
	A = R G
	F = A^3
	e0 = F_0
	e1 = F_1
	a = e0*e*f
	b = e1*e*f
	c = e0*e*f*g	  
        d=a+b+c
	ac=allComponent(d)

///

doc ///
   Key
     pathTerms
     (pathTerms,PAElement)
   Headline
     Get pairs of coefficient and monomial of a Path Algebra Element
   Usage
     pt=pathTerms(f)
   Inputs
     f:PAElement
   Outputs
     pt:List
   Description
      Text 
       This function returns a List of pairs of coefficient and monomial of a Path Algebra Element
      Example
        adj = matrix {{1,0},{1,1}}
	G = paGraph({v,w},{e,f,g},adj)
	R = QQ
	A = R G
	a = 2*e*f+ 3*e*f*g	  
	pt=pathTerms(a)

///

doc ///
   Key
     (pathTerms,PAVector)
   Headline
     Get pairs of coefficient and monomial of a Path Algebra Vector
   Usage
     pt=pathTerms(f)
   Inputs
     f:PAVector
   Outputs
     pt:List
   Description
      Text 
       This function returns a List of pairs of coefficient and module monomial of a Path Algebra Vector
      Example
        adj = matrix {{1,0},{1,1}}
	G = paGraph({v,w},{e,f,g},adj)
	R = QQ
	A = R G
	F = A^3
	e0 = F_0
	e1 = F_1
	a = e0*e*f
	b = e1*e*f
	c = e0*e*f*g	  
        d=a+b+c
        pt=pathTerms(d)

///


doc ///
   Key
     selectComponent
     (selectComponent,PAVector,ZZ)
   Headline
     Get terms of a PAVector of a certain component
   Usage
     sc=selectComponent(f,x)
   Inputs
     f:PAVector
     x:ZZ
   Outputs
     sc:PAVector
   Description
      Text 
       This function returns a PAVector of Path Algebra Modulde Monomial of a certain component
      Example
        adj = matrix {{1,0},{1,1}}
	G = paGraph({v,w},{e,f,g},adj)
	R = QQ
	A = R G
	F = A^3
	e0 = F_0
	e1 = F_1
	a = e0*e*f
	b = e1*e*f
	c = e0*e*f*g	  
        d=a+b+c
	selectComponent(d,0)
	selectComponent(d,1)

///

doc ///
   Key
      leadPair
      (leadPair, PAVector)
   Headline
      Returns the lead term and the coefficient of an PAVector
   Usage
     lp = leadPair f
   Inputs
     f : PAVector
   Outputs
     lp : Sequence
   Description
      Text
        Returns the lead term and the coefficient of an PAVector
      Example
        adj = matrix {{1,0},{1,1}}
	G = paGraph({v,w},{e,f,g},adj)
	R = QQ
	A = R G
	F = A^3
	e0 = F_0
	e1 = F_1
	a = e0*e*f
	b = e1*e*f
	c = e0*e*f*g	  
        d=a+b+c
	lp=leadPair(d)
///

doc ///
   Key
      isVertex
      (isVertex, PAPath)
   Headline
      Determine whether the input is a vertex
   Usage
     iv = isVertex f
   Inputs
     f : PAPath
   Outputs
     iv : Boolean
   Description
      Text
        Determine whether a Path is a vertex
      Example
        adj = matrix {{1,0},{1,1}}
	G = paGraph({v,w},{e,f,g},adj)
	a=paPath(G,{0})
	isVertex(a)
	b=paPath(G,{})
	isVertex(b)

///

doc ///
   Key
      (isVertex, PAElement)
   Headline
      Determine whether the input is a vertex
   Usage
     iv = isVertex f
   Inputs
     f : PAElement
   Outputs
     iv : Boolean
   Description
      Text
        Determine whether a Path Algebra Element is a vertex
      Example
        adj = matrix {{1,0},{1,1}}
	G = paGraph({v,w},{e,f,g},adj)
	R = QQ
	A = R G
        a=putInPathAlgebra(A,paPath(G,{0,0,1,2}))	
	isVertex(a)
	b=w
	isVertex(b)
	c=v
	isVertex(c)
	d=e*f*g
	isVertex(d)

///

doc ///
   Key
      (isVertex, PAVector)
   Headline
      Determine whether the input is a vertex
   Usage
     iv = isVertex f
   Inputs
     f : PAVector
   Outputs
     iv : Boolean
   Description
      Text
        Determine whether a Path Algebra Vector is a vertex
      Example
        adj = matrix {{1,0},{1,1}}
	G = paGraph({v,w},{e,f,g},adj)
	R = QQ
	A = R G
	F = A^3
	a = (F_0)*e	
	isVertex(a)
	b = (F_0)*v
	isVertex(b)

///

doc ///
   Key
      putInPathAlgebra
      (putInPathAlgebra,PathAlgebra,PAPath)
   Headline
      Create an Path Algebra
   Usage
      M = putInPathAlgebra(K,f)
   Inputs
      K : PathAlgebra
      f : PAPath
   Outputs
      M : PathAlgebra
   Description
      Text
         This command put a PAPath into a PathAlgebra.
      Example
         M = matrix {{1,0},{1,1}}
	 G = paGraph({v,w},{e,f,g},M,Weights=>{2,1,1},Degrees => {2,1,1})
	 R = QQ
	 A = R G
	 putInPathAlgebra(A,paPath(G,{0,1}))
///

      
doc ///
   Key
      (putInPathAlgebra,PathAlgebra,PAPath,QQ)
   Headline
      Create an Path Algebra
   Usage
      M = putInPathAlgebra(K,f,x)
   Inputs
      K : PathAlgebra
      f : PAPath
      x : QQ
   Outputs
      M : PathAlgebra
   Description
      Text
         This command put a PAPath into a PathAlgebra.
      Example
         M = matrix {{1,0},{1,1}}
	 G = paGraph({v,w},{e,f,g},M,Weights=>{2,1,1},Degrees => {2,1,1})
	 R = QQ
	 A = R G
	 putInPathAlgebra(A,paPath(G,{0,1}),3_R)
///

doc ///
   Key
      (putInPathAlgebra,PathAlgebra,PAPath,ZZ)
   Headline
      Create an Path Algebra
   Usage
      M = putInPathAlgebra(K,f,x)
   Inputs
      K : PathAlgebra
      f : PAPath
      x : ZZ
   Outputs
      M : PathAlgebra
   Description
      Text
         This command put a PAPath into a PathAlgebra.
      Example
         M = matrix {{1,0},{1,1}}
	 G = paGraph({v,w},{e,f,g},M,Weights=>{2,1,1},Degrees => {2,1,1})
	 R = ZZ
	 A = R G
	 putInPathAlgebra(A,paPath(G,{0,1}),3_ZZ)
///


doc ///
   Key
      (putInPathAlgebra,PathAlgebra,PAVector)
   Headline
      Create an Path Algebra
   Usage
      M = putInPathAlgebra(A,f)
   Inputs
      A : PathAlgebra
      f : PAVector
   Outputs
      M : PathAlgebra
   Description
      Text
         This command put a PAVector into a PathAlgebra.
      Example	 
     	 adj = matrix {{3}}
	 G = paGraph({v},{x,y,z},adj)
	 R = QQ
	 A = R G
	 F0 = A^1
	 e0 = F0_0
	 f = e0*x*y*z
	 putInPathAlgebra(A,f)	 

///



doc ///
   Key
      (putInPathAlgebra,PathAlgebra,List)
   Headline
      Create a list of Path Algebra
   Usage
      M = putInPathAlgebra(A,L)
   Inputs
      A : PathAlgebra
      L : List
   Outputs
      M : List 
   Description
      Text
         This command put a list of PAVector into a list of PathAlgebra.
      Example
     	 adj = matrix {{3}}
	 G = paGraph({v},{x,y,z},adj)
	 R = QQ
	 A = R G
	 F0 = A^1
	 e0 = F0_0
	 f = e0*x*y*z
	 g = e0*x*y*y
	 h = e0*y*z
	 L = {f,g,h}
	 putInPathAlgebra(A,L)	
///

doc ///
   Key
      (putInPathAlgebra,PathAlgebra,PAPath,RingElement)
   Headline
      Create an Path Algebra
   Usage
      M = putInPathAlgebra(A,f,r)
   Inputs
      A : PathAlgebra
      f : PAPath
      r : RingElement
   Outputs
      M : PathAlgebra
   Description
      Text
         This command put a PAVector into a PathAlgebra.
      Example	 
     	 adj = matrix {{3}}
	 G = paGraph({v},{x,y,z},adj)
	 R = QQ
	 A = R G
	 r = 1_R
	 f = paPath(G,{0,1,2})
	 putInPathAlgebra(A,f,r)	 

///

doc ///
   Key
      (putInPathAlgebra,PathAlgebraQuotient,PAPath)
   Headline
      Create a list of Path Algebra Quotient
   Usage
      M = putInPathAlgebra(A,L)
   Inputs
      A : PathAlgebraQuotient
      L : PAPath
   Outputs
      M : PathAlgebra 
   Description
      Text
         This command put a PAPath into a PathAlgebraQuotient.
      Example
         M = matrix {{3}}
	 G = paGraph({v},{a,b,c},M,Weights=>{1,1,1})
	 kk = ZZ/32003
	 R = kk G
	 I = paIdeal {2*a*b + 3*b*a + 5*c^2,
                      2*b*c + 3*c*b + 5*a^2,
                      2*c*a + 3*a*c + 5*b^2}
	 S = R/I
	 x = a^2
	 y = putInPathAlgebra(S,paPath(G,{0,0}))
         x == y
///


doc ///
   Key
      (putInPathAlgebra,PathAlgebraQuotient,PAPath,RingElement)
   Headline
      Create an Path Algebra
   Usage
      M = putInPathAlgebra(K,f,x)
   Inputs
      K : PathAlgebraQuotient
      f : PAPath
      x : RingElement
   Outputs
      M : PathAlgebraQuotient
   Description
      Text
         This command put a PAPath into a PathAlgebraQuotient.
      Example
         M = matrix {{3}}
	 G = paGraph({v},{a,b,c},M,Weights=>{1,1,1})
	 kk = ZZ/32003
	 R = kk G
	 I = paIdeal {2*a*b + 3*b*a + 5*c^2,
                      2*b*c + 3*c*b + 5*a^2,
                      2*c*a + 3*a*c + 5*b^2}
	 S = R/I
	 y = putInPathAlgebra(S,paPath(G,{0,0}),3_(kk))
///

doc ///
   Key
      isSubModMon
      (isSubModMon,PAElement,PAModMon)
   Headline
      A prefix check for Path Algebra Module Monomial
   Usage
      C = isSubModMon(f,g)
   Inputs
      f : PAElement
      g : PAModMon
   Outputs
      C : Boolean
   Description
      Text
         This command checks if a PAElement is a prefix of a Path Algebra Module Monomial
      Example
	 adj = matrix {{1,0},{1,1}}
	 G = paGraph({v,w},{e,f,g},adj)
	 R = QQ
	 A = R G
         L = e*f
	 H = paModMon(0,paPath(G,{0,0,1,2}))
	 K = e*e
	 isSubModMon(L,H)
	 isSubModMon(K,H)

///

doc ///
   Key
      (isSubModMon,PAModMon,PAModMon)
   Headline
      A prefix check for Path Algebra Module Monomial
   Usage
      C = isSubModMon(f,g)
   Inputs
      f : PAModMon
      g : PAModMon
   Outputs
      C : Boolean
   Description
      Text
         This command checks if a Path Algebra Module Monomial is a prefix of another
      Example
	 adj = matrix {{1,0},{1,1}}
	 G = paGraph({v,w},{e,f,g},adj)
	 R = QQ
	 A = R G
	 L = paModMon(0,paPath(G,{0,1}))
	 H = paModMon(0,paPath(G,{0,0,1,2}))
	 K = paModMon(0,paPath(G,{0,0}))
	 isSubModMon(L,H)
	 isSubModMon(K,H)

///

doc ///
   Key
      (isSubModMon,PAModMon,PAPath)
   Headline
      A prefix check for Path Algebra Module Monomial
   Usage
      C = isSubModMon(f,g)
   Inputs
      f : PAModMon
      g : PAPath
   Outputs
      C : Boolean
   Description
      Text
         This command checks if a Path Algebra Module Monomial is a prefix of a PAPath
      Example
	 adj = matrix {{1,0},{1,1}}
	 G = paGraph({v,w},{e,f,g},adj)
	 R = QQ
	 A = R G
	 L = paModMon(0,paPath(G,{0,1}))
	 H = paPath(G,{0,0,1,2})
	 K = paModMon(0,paPath(G,{0,0}))
	 isSubModMon(L,H)
	 isSubModMon(K,H)

///


doc ///
   Key
      (isSubModMon,PAModMon,PAElement)
   Headline
      A prefix check for Path Algebra Module Monomial
   Usage
      C = isSubModMon(f,g)
   Inputs
      f : PAModMon
      g : PAElement
   Outputs
      M : Boolean
   Description
      Text
         This command checks if a Path Algebra Module Monomial is a prefix of a Path Algebra Element
      Example
	 adj = matrix {{1,0},{1,1}}
	 G = paGraph({v,w},{e,f,g},adj)
	 R = QQ
	 A = R G
	 L = paModMon(0,paPath(G,{0,1}))
	 H = e*e*f*g
	 K = paModMon(0,paPath(G,{0,0}))
	 isSubModMon(L,H)
	 isSubModMon(K,H)

///

doc ///
   Key
      allSubModMon
      (allSubModMon,PAVector,List)	      
   Headline
      A prefix check for Path Algebra Vector
   Usage
      C = allSubModMon(f,L)
   Inputs
      f : PAVector
      L : List
   Outputs
      C : List
   Description
      Text
         "To DO" 
      Example
	 adj = matrix {{1,0},{1,1}}
	 G = paGraph({v,w},{e,f,g},adj)
	 R = QQ
	 A = R G


///


doc ///
   Key
      stdBasisVector
      (stdBasisVector,PAModule,ZZ)
   Headline
      Find the stand basis vector of a Path Algebra Module
   Usage
      V = stdBasisVector(M,n)
   Inputs
      M : PAModule
      n : ZZ
   Outputs
      M : PAVector 
   Description
      Text
         Find the stand basis vector of a Path Algebra Module with a given module index
      Example
     	 adj = matrix {{3}}
	 G = paGraph({v},{x,y,z},adj)
	 R = QQ
	 A = R G
	 F3 = A^3
	 stdBasisVector(F3,2)
	
///

doc ///
   Key
      removeZeropath 
      (removeZeropath,List)
   Headline
      Remove all of the zero paths
   Usage
      R = removeZeropath(L)
   Inputs
      L : List
   Outputs
      R : List
///


doc ///
   Key
      makeMonic
      (makeMonic,PAElement)
   Headline
      Make a Path Algebra Element Monic
   Usage
      mm = makeMonic(f)
   Inputs
      f : PAElement
   Outputs
      mm : PAElement
   Description
      Text
        Make a Path Algebra Element Monic
      Example
     	 adj = matrix {{3}}
	 G = paGraph({v},{x,y,z},adj)
	 R = QQ
	 A = R G
	 f = 5*x*y*z+3*y*z
	 g = makeMonic(f)
	
///

doc ///
   Key
      (makeMonic,PAVector)
   Headline
      Make a PAVector Monic
   Usage
      mm = makeMonic(f)
   Inputs
      f : PAVector
   Outputs
      mm : PAVector
   Description
      Text
         Make a PAVector Monic
      Example
     	 adj = matrix {{3}}
	 G = paGraph({v},{x,y,z},adj)
	 R = QQ
	 A = R G
	 F = A^3
	 f = (F_0)*5*x*y*z+(F_1)*3*y*z
	 g = makeMonic(f)
	
///
 	     

doc ///
   Key
      isPrefix
      (isPrefix,List,List)
   Headline
      A prefix check for list of numbers
   Usage
      ispre = isPrefix(L,M)
   Inputs
      L : List
      M : List
          L and M are two lists of numbers;
   Outputs
      ispre : Boolean
   Description
      Text
        Check if the first list of numbers is a prefix of the second list.
      Example
        L = {0,1,2,3,4}
	M = {0,1,2,3,4,4,3,2,1}
	isPrefix(L,M)
	isPrefix(M,L)
	
///

doc ///
   Key
      (isPrefix,PAModMon,PAModMon)
   Headline
      A prefix check for Path Algebra Module Monomials 
   Usage
      ispre = isPrefix(f,g)
   Inputs
      f : PAModMon
      g : PAModMon

   Outputs
      ispre : Boolean
   Description
      Text
        Check if the first PAModMon is a prefix of the second.
      Example
       "TODO"
	
///

doc ///
   Key
      (isPrefix,PAVector,PAVector)
   Headline
      A prefix check for PAVectors
   Usage
      ispre = isPrefix(f,g)
   Inputs
      f : PAVector
      g : PAVector

   Outputs
      ispre : Boolean
   Description
      Text
        Check if the first PAVector is a prefix of the second.
      Example
       "TODO"
	
///


doc ///
   Key
      findSubModMon
      (findSubModMon,PAVector,List,List)
   Headline
      A subModMon check for PAVectors
   Usage
      fsmm = findSubModMon(f,L,M)
   Inputs
      f : PAVector
      L : List
      M : List
   Outputs
      fsmm : Sequence
   Description
      Text
        Check if the PAVector is a prefix of either list of elements.
        The first list are module generators and the second list are defining equations of the algebra.
	The function checks the defining equations of the algebras first.
      Example
        "TO DO"
	
///

doc ///
   Key
      (findSubModMon,Sequence,List,List)
   Headline
      A subModMon check for PAVectors
   Usage
      fsmm = findSubModMon(f,L,M)
   Inputs
      f : Sequence
      L : List
      M : List
   Outputs
      fsmm : Sequence
   Description
      Text
        f is a tuple which the first entry is a PAVector and the second entry is track info.
        This function checks if the PAVector is a prefix of the list of module generators, or if any of the list of defining equations of the algebras is a suffix of the PAVector
	The function checks the defining equations of the algebras first.
      Example
        "TO DO"	
///


doc ///
   Key
      divAlgorithm
      (divAlgorithm,List,PAElement)
   Headline
      A division algorithm for PAElements
   Usage
      divAl = divAlgorithm(L,f)
   Inputs
      f : PAElement
      L : List
   Outputs
      divAl : Sequence
   Description
      Text
        This functions returns a remainder of f divided by a list of PAPaths.
	The first entry of return value gives the remainder.
      Example
        "TO DO"
	
///

doc ///
   Key
      numEdges
      (numEdges,PAGraph)
   Headline
      Find the number of Edges of a PAGraph
   Usage
      numE = numEdges(G)
   Inputs
      G : PAGraph
   Outputs
      numE : ZZ
   Description
      Text
        This functions returns the number of Edges of a PAGraph
      Example
        M = matrix {{1,0},{1,1}}
	G = paGraph({v,w},{e,f,g},M)
	numEdges(G)
	
///

doc ///
   Key
      (numEdges,PathAlgebra) 
   Headline
      Find the number of Edges of a Path Algebra
   Usage
      numE = numEdges(A)
   Inputs
      A : PathAlgebra
   Outputs
      numE : ZZ
   Description
      Text
        This functions returns the number of Edges of a Path Algebra
      Example
        M = matrix {{1,0},{1,1}}
	G = paGraph({v,w},{e,f,g},M)
	R = QQ
	A = R G
	numEdges(A)
	
///

doc ///
   Key
      (numEdges,PathAlgebraQuotient) 
   Headline
      Find the number of Edges of a Path Algebra Quotient
   Usage
      numE = numEdges(A)
   Inputs
      A : PathAlgebraQuotient
   Outputs
      numE : ZZ
   Description
      Text
        This functions returns the number of Edges of a Path Algebra Quotient
      Example
        M = matrix {{1,0},{1,1}}
	G = paGraph({v,w},{e,f,g},M)
	R = QQ
	A = R G
	I=paIdeal{e*f,f*g}
	B = A/I
	numEdges(B)
	
///

doc ///
   Key
      numVertices
      (numVertices,PAGraph)
   Headline
      Find the number of Vertices of a PAGraph
   Usage
      numV =  numVertices(G)
   Inputs
      G : PAGraph
   Outputs
      numE : ZZ
   Description
      Text
        This functions returns the number of Vertices of a PAGraph
      Example
        M = matrix {{1,0},{1,1}}
	G = paGraph({v,w},{e,f,g},M)
	numVertices(G)
	
///

doc ///
   Key
      (numVertices,PathAlgebra) 
   Headline
      Find the number of Vertices of a Path Algebra
   Usage
      numE = numVertices(A)
   Inputs
      A : PathAlgebra
   Outputs
      numE : ZZ
   Description
      Text
        This functions returns the number of Vertices of a Path Algebra
      Example
        M = matrix {{1,0},{1,1}}
	G = paGraph({v,w},{e,f,g},M)
	R = QQ
	A = R G
	numVertices(A)
	
///

doc ///
   Key
      (numVertices,PathAlgebraQuotient) 
   Headline
      Find the number of Vertices of a Path Algebra Quotient
   Usage
      numE = numVertices(A)
   Inputs
      A : PathAlgebraQuotient
   Outputs
      numE : ZZ
   Description
      Text
        This functions returns the number of Vertices of a Path Algebra Quotient
      Example
        M = matrix {{1,0},{1,1}}
	G = paGraph({v,w},{e,f,g},M)
	R = QQ
	A = R G
	I=paIdeal{e*f,f*g}
	B = A/I
	numVertices(B)
	
///

doc ///
  Key
     (use, PathAlgebra)
     (use, PathAlgebraQuotient)
  Headline
     Brings the variables of a particular PathAlgebra in scope
  Usage
    use A
  Inputs
    A : PathAlgebra
  Description
    Text
       This function brings the variables of a particular PathAlgebra in scope.
       For an illustration:
    Example
       adj = matrix {{3}}
       G = paGraph({v},{x,y,z},adj)
       R = QQ
       A = R G
       x
       I = {x*y - y*x, x*z - z*x, y*z - z*y}
       I = paIdeal I
       buchAlgorithm I
       B = A/I
       x
    Text
       As you can see, at this point the interpreter treats x elements of B.  To go back to
       A, we run the command use A:
    Example
       use A
       x
    Text
       Similarly, to go back to B, we run the command use B:
    Example
       use B
       x
///


doc ///
   Key
      (isWellDefined, PAMap)
   Headline
      Determines if an PAMap is well-defined.
   Usage
      isWellDefined f
   Inputs
      f : PAMap
   Outputs
      : Boolean
   Description
      Text
         Returns true if the given PAMap evaluates as 0 on the defining relations
	 of the source.
      Example
         adj = matrix {{3}}
	 G = paGraph({v},{x,y,z},adj)
	 R = QQ
	 A = R G
	 F0 = A^1
	 I = {x*y - y*x, x*z - z*x, y*z - z*y}
	 I = paIdeal I
	 buchAlgorithm I
	 B = A/I
	 phi = paMap(B,A,{v},{x,y,z})
	 isWellDefined phi
	 use A
	 psi = paMap(A,B,{v},{x,y,z})
	 isWellDefined psi

///

doc ///
   Key
     (ambient,PAMap)
   Headline
     Ambient ring of an PAMap
   Usage
     A = ambient B 
   Inputs
     B : PAMap
   Outputs
     A : PathAlgebra
   Description
      Text
         Returns the ambient ring of an @ TO PAMap  @.  
	 
	 As quotients of PathAlgebras are added, this will return the top-level ambient ring.
	 
      Example
         "TODO"
///

doc ///
   Key
     (ambient,PathAlgebraQuotient)
   Headline
     Ambient ring of an PathAlgebraQuotient
   Usage
     A = ambient B 
   Inputs
     B : PathAlgebraQuotient
   Outputs
     A : PathAlgebra
   Description
      Text
         Returns the ambient ring of an @ TO PathAlgebraQuotient  @.  
	 
	 As quotients of PathAlgebras are added, this will return the top-level ambient ring.
	 
      Example
         "TODO"
///

doc ///
  Key
    buchPAModule
    (buchPAModule,List,List)
  Headline
    Computes the Grobner Basis of a list of PAVectors
  Usage
    G=buchPAModule(M,I)
  Inputs
    M:List
    I:List
  Outputs
    G:List
  Description
    Text
     "to do"
    Example
     "TODO"
///
