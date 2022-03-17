/**
@file parser.hh
@brief Functions to import Grid from grid files and export grid structures and vectors to vtk files
*/

#ifndef PARSER_HH
#define PARSER_HH

namespace conformingsimplexgrid {

    template< class Grid >
    void importFromFile( Grid &T, const char* file )
    /// import elements and coordinates from grid file.
    {
        static const int NV = Grid::NV;
        static const int dim = Grid::dim;
        FILE *fp;
        int numElements, numNodes;
        char puffer[30];

        fp = fopen(file, "r");
        if(fp == NULL) {
            printf("Fehler bei fopen() von %s... \n", file);
            abort();
        }
        
        if (fscanf(fp, "nE = %d; nC = %d;\n", &numElements, &numNodes))
        {
        
            T.resizeElements( numElements );
            T.resizeNodes( numNodes );

            rewind( fp );
            unsigned int temp = 0;
            while( fgets(puffer, 10, fp) != NULL ) {
                if(strstr(puffer,"n4e") != 0){
                    for(size_t i = 0; i < T.nElements(); i++)
                        for(size_t j = 0; j < NV; j++){
                            if (!fscanf(fp,"%d",&temp)) {printf("Fehler bei fscanf()\n");}
                            T.elements(i,j) = temp; 
                        }
                }
            }
            rewind( fp );
            double temp_double = 0.;
            while( fgets(puffer, 10, fp) != NULL ) {
                if(strstr(puffer,"c4n") != 0){
                    for(size_t i=0; i < T.nNodes(); i++)
                        for(size_t j=0; j < dim; j++){
                            if (!fscanf(fp,"%lf",&temp_double)) {printf("Fehler bei fscanf()\n");}
                            T.coordinates(i,j) = temp_double;
                        }
                }
            }
            
            T.initialTriangulation();
        }
        else
        {
            printf("Fehler bei fscanf()\n");
            abort();
        }
    }

    template< class Grid, class Array>
    void export2vtk( Grid & T, const Array & u, const char *file )
    ///export grid structure and solution vector u to a vtk file.
    {
        FILE *fid;
        fid = fopen(file,"w");  
        const size_t NV = Grid::NV;
        const size_t dimworld = Grid::dim;
        if( fid == NULL ) {
            printf( "Fehler bei fopen() von '%s'\n", file );
            abort();
        }
        fprintf(fid,"<?xml version=\"1.0\"?>\n");
        fprintf(fid,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
        fprintf(fid,"  <UnstructuredGrid>\n");
        fprintf(fid,"    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", (int)T.nNodes() ,(int)T.nElements() );
        fprintf(fid,"      <PointData Scalars=\"scalars\">\n");      
        fprintf(fid,"        <DataArray type=\"Float32\" Name=\"phase\" format=\"ascii\">\n");  
        for (size_t i = 0; i < T.nNodes();i++)
            fprintf(fid,"          %3.3f\n", u[i]);     
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"      </PointData>\n");
        fprintf(fid,"      <Points>\n");
        fprintf(fid,"        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for (size_t i = 0; i < T.nNodes(); i++){
            fprintf(fid,"         ");
            for (size_t j = 0; j < dimworld; j++)
                fprintf(fid," %3.9f", T.coordinates(i,j));
            if(dimworld == 2)
                fprintf(fid," %3.9f", u[i]);
            fprintf(fid,"\n");  
        }
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"      </Points>\n");
        fprintf(fid,"      <Cells>\n");
        fprintf(fid,"        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
        for (size_t i = 0; i < (size_t)T.nElements(); i++){
            fprintf(fid,"         ");
            for (size_t j = 0; j < NV; j++)
                fprintf(fid," %d", (int)T.elements(i,j));
            fprintf(fid,"\n");  
        }       
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
        fprintf(fid,"          ");      
        for (size_t i = NV; i < NV*T.nElements()+1; i+=NV)
                fprintf(fid,"%d\n",(int)i);
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
        for (size_t i=0; i < T.nElements();i++)
                        fprintf(fid,"          %d\n",(int)((NV-1)*(NV-1)+1));
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"      </Cells>\n");
        fprintf(fid,"    </Piece>\n");
        fprintf(fid,"  </UnstructuredGrid>\n");
        fprintf(fid,"</VTKFile>\n");
        fclose(fid);
    }
    
    
    template< class Grid, class Array, class ElementVector>
    void export2vtk( Grid & T, const Array & u, ElementVector & partition, ElementVector & markedElements, const char *file )
    ///export grid structure, solution vector u, partition and markedElements vector to a vtk file.
    {
        FILE *fid;
        fid = fopen(file,"w");	
        const size_t NV = Grid::NV;
        const size_t dimworld = Grid::dim;
        if( fid == NULL ) {
            printf( "Fehler bei fopen() von '%s'\n", file );
            abort();
        }
        fprintf(fid,"<?xml version=\"1.0\"?>\n");
        fprintf(fid,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
        fprintf(fid,"  <UnstructuredGrid>\n");
        fprintf(fid,"    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", (int)T.nNodes() ,(int)T.nElements() );
        fprintf(fid,"      <PointData Scalars=\"scalars\">\n");      
        fprintf(fid,"        <DataArray type=\"Float32\" Name=\"phase\" format=\"ascii\">\n");	
        for (size_t i = 0; i < T.nNodes();i++)
            fprintf(fid,"          %3.3f\n", u[i]);	
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"      </PointData>\n");
        fprintf(fid,"      <CellData Scalars=\"scalars\">\n");
        fprintf(fid,"        <DataArray type=\"UInt8\" Name=\"partition\" format=\"ascii\">\n");
        for (size_t i=0; i < T.nElements();i++)
           fprintf(fid," %d\n",partition[i]);
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"      </CellData>\n");
        fprintf(fid,"      <CellData Scalars=\"scalars\">\n");
        fprintf(fid,"        <DataArray type=\"UInt8\" Name=\"markedElements\" format=\"ascii\">\n");
        for (size_t i=0; i < T.nElements();i++)
           fprintf(fid," %d\n",markedElements[i]);
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"      </CellData>\n");   
        fprintf(fid,"      <Points>\n");
        fprintf(fid,"        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for (size_t i = 0; i < T.nNodes(); i++){
            fprintf(fid,"         ");
            for (size_t j = 0; j < dimworld; j++)
                fprintf(fid," %3.9f", T.coordinates(i,j));
            if(dimworld == 2)
                fprintf(fid," %3.9f", u[i]);
            fprintf(fid,"\n");	
        }
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"      </Points>\n");
        fprintf(fid,"      <Cells>\n");
        fprintf(fid,"        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
        for (size_t i = 0; i < (size_t)T.nElements(); i++){
            fprintf(fid,"         ");
            for (size_t j = 0; j < NV; j++)
                fprintf(fid," %d", (int)T.elements(i,j));
            fprintf(fid,"\n");	
        }	
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
        fprintf(fid,"          ");	
        for (size_t i = NV; i < NV*T.nElements()+1; i+=NV)
                fprintf(fid,"%d\n",(int)i);
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
        for (size_t i=0; i < T.nElements();i++)
                        fprintf(fid,"          %d\n",(int)((NV-1)*(NV-1)+1));
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"      </Cells>\n");
        fprintf(fid,"    </Piece>\n");
        fprintf(fid,"  </UnstructuredGrid>\n");
        fprintf(fid,"</VTKFile>\n");
        fclose(fid);
    }
    
    template< class Grid, class Array, class ElementVector>
    void export2vtk_faces( Grid & T, const Array & u, ElementVector & partition, ElementVector & markedElements, const char *file )
    ///export the faces grid, solution vector u, partition and markedElements vector to a vtk file.
    {
        FILE *fid;
        fid = fopen(file,"w");  
        const size_t NV = Grid::NV-1;
        const size_t dimworld = Grid::dim;
        
        if( fid == NULL ) {
            printf( "Fehler bei fopen() von '%s'\n", file );
            abort();
        }

        fprintf(fid,"<?xml version=\"1.0\"?>\n");
        fprintf(fid,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
        fprintf(fid,"  <UnstructuredGrid>\n");
        fprintf(fid,"    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", (int)T.nNodes() ,(int)T.nFaces() );
        fprintf(fid,"      <PointData Scalars=\"scalars\">\n");      
        fprintf(fid,"        <DataArray type=\"Float32\" Name=\"phase\" format=\"ascii\">\n");  
        for (size_t i = 0; i < T.nNodes();i++)
            fprintf(fid,"          %3.3f\n", u[i]);     
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"      </PointData>\n");
        fprintf(fid,"      <CellData Scalars=\"scalars\">\n");
        fprintf(fid,"        <DataArray type=\"UInt8\" Name=\"partition\" format=\"ascii\">\n");
        for (size_t i=0; i < T.nFaces();i++)
           fprintf(fid," %d\n",partition[i]);
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"      </CellData>\n");
        fprintf(fid,"      <CellData Scalars=\"scalars\">\n");
        fprintf(fid,"        <DataArray type=\"UInt8\" Name=\"markedElements\" format=\"ascii\">\n");
        for (size_t i=0; i < T.nFaces();i++)
           fprintf(fid," %d\n",markedElements[i]);
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"      </CellData>\n");   
        fprintf(fid,"      <Points>\n");
        fprintf(fid,"        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for (size_t i = 0; i < T.nNodes(); i++){
            fprintf(fid,"         ");
            for (size_t j = 0; j < dimworld; j++)
                fprintf(fid," %3.9f", T.coordinates(i,j));
            if(dimworld == 2)
                fprintf(fid," %3.9f", 0.0);
            fprintf(fid,"\n");  
        }
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"      </Points>\n");
        fprintf(fid,"      <Cells>\n");
        fprintf(fid,"        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
        for (size_t i = 0; i < (size_t)T.nFaces(); i++){
            fprintf(fid,"         ");
            for (size_t j = 0; j < NV; j++)
                fprintf(fid," %d", (int)T.faces(i,j));
            fprintf(fid,"\n");  
        }       
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
        fprintf(fid,"          ");      
        for (size_t i = NV; i < NV*T.nFaces()+1; i+=NV)
                fprintf(fid,"%d\n",(int)i);
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
        for (size_t i=0; i < T.nFaces();i++)
                        fprintf(fid,"          %d\n",(int)((NV-1)*(NV-1)+1));
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"      </Cells>\n");
        fprintf(fid,"    </Piece>\n");
        fprintf(fid,"  </UnstructuredGrid>\n");
        fprintf(fid,"</VTKFile>\n");
        fclose(fid);
    }
    
#if 0
   template< int NV, int dimworld >
   void export2vtu( const mesh< NV, dimworld > *T, const double *u, const double *part, const char *file ){
        int i, j;
        FILE *fid;
        fid = fopen(file,"w");  

        if( fid == NULL ) {
            printf( "Fehler bei fopen() von '%s'\n", file );
            abort();
        }

        fprintf(fid,"<?xml version=\"1.0\"?>\n");
        fprintf(fid,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
        fprintf(fid,"  <UnstructuredGrid>\n");
        fprintf(fid,"    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",(int)T->c4n.size(),(int)T->n4e.size());
        fprintf(fid,"      <PointData Scalars=\"scalars\">\n");      
        fprintf(fid,"        <DataArray type=\"Float32\" Name=\"phase\" format=\"ascii\">\n");  
        for (i=0;i<(int)T->c4n.size();i++)
            fprintf(fid,"          %3.3f\n", u[i]); 
        fprintf(fid,"        </DataArray>\n");     
        fprintf(fid,"        <DataArray type=\"Float32\" Name=\"partition\" format=\"ascii\">\n");  
        for (i=0;i<(int)T->c4n.size();i++)
            fprintf(fid,"          %3.3f\n", part[i] + 1.); 
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"      </PointData>\n");
        fprintf(fid,"      <Points>\n");
        fprintf(fid,"        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for (i=0;i<(int)T->c4n.size();i++){
            fprintf(fid,"         ");
            for (j=0;j<dimworld;j++)
                fprintf(fid," %3.9f", T->c4n(i,j));
            if(dimworld == 2)
                fprintf(fid," %3.9f", 0.0);
            fprintf(fid,"\n");  
        }
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"      </Points>\n");
        fprintf(fid,"      <Cells>\n");
        fprintf(fid,"        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
        for (i=0;i<(int)T->n4e.size();i++){
            fprintf(fid,"         ");
            for (j=0;j<NV;j++)
                fprintf(fid," %d", (int)T->n4e(i,j));
            fprintf(fid,"\n");  
        }   
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
        fprintf(fid,"          ");  
        for (i=NV;i<(NV)*(int)T->n4e.size()+1; i+=NV)
            fprintf(fid,"%d\n",i);
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
        for (i=0;i<(int)T->n4e.size();i++)
                fprintf(fid,"          %d\n",(NV-1)*(NV-1)+1);
        fprintf(fid,"        </DataArray>\n");
        fprintf(fid,"      </Cells>\n");
        fprintf(fid,"    </Piece>\n");
        fprintf(fid,"  </UnstructuredGrid>\n");
        fprintf(fid,"</VTKFile>\n");
        fclose(fid);
    }


    template< int NV, int dimworld >
    void export2matlab( mesh<NV,dimworld> &T, ManagedArray<double, NV*dimworld> normals, char *file ){
	    int i, k;
	    FILE *fid;
	    fid = fopen(file,"w");	
	    fprintf(fid,"c4n = [ ");
	    for (i=0;i<(int)T.c4n.size();i++){
		    for(k=0;k<dimworld;k++)
			    fprintf(fid," %3.3f ", T.c4n(i,k));
		    fprintf(fid,";...\n");
	    } fprintf(fid,"];\n");
	    fprintf(fid,"n4e = [ ");
	    for (i=0;i<(int)T.n4e.size();i++){
		    for(k=0;k<NV;k++)
			    fprintf(fid," %d ", T.n4e(i,k));
		    fprintf(fid,";...\n");
	    } fprintf(fid,"];\n");
	    fprintf(fid,"normals = [ ");
	    for (i=0;i<(int)T.n4e.size();i++){
		    for(k=0;k<dimworld*NV;k++)
			    fprintf(fid," %3.3f ", normals(i,k));
		    fprintf(fid,";...\n");
	    } fprintf(fid,"];\n");

	    fclose(fid);
    }

#endif 
} // namespace conformingsimplexgrid

#endif // PARSER_HH