//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Utilities;

import com.Sts.Actions.Wizards.VirtualVolume.*;

public class StsConvolve
{
    public static final String BARTLETT = StsFilterVVolumePanel.BARTLETT;
    public static final String GAUSSIAN = StsFilterVVolumePanel.GAUSSIAN;
    public static final String LAPLACIAN = StsFilterVVolumePanel.LAPLACIAN;

    public static final int XDIR = StsParameters.XDIR;
    public static final int YDIR = StsParameters.YDIR;
    public static final int ZDIR = StsParameters.ZDIR;

   static public float[] convolveFloat3D(int dir, int nPlaneRows, int nPlaneCols, int nPlaneMin, int nPlaneMax,
                                         float[][] inPlaneData, String kernelName, int xSize, int ySize, int zSize)
    {
       Kernel kernel = new Kernel(dir, kernelName, xSize, ySize, zSize);

       int nPlanes = nPlaneMax - nPlaneMin + 1;
       int rowSize = getRowSize(dir, ySize, zSize);
       int colSize = getColSize(dir, xSize, ySize);

       return convolveFloat3D(inPlaneData, nPlanes, nPlaneRows, nPlaneCols, rowSize, colSize, kernel);
   }

    static public float[] convolveFloat3D(float[][]inPlaneData, int nPlanes, int nPlaneRows, int nPlaneCols,
                                          int rowSize, int colSize, Kernel kernel)
    {
        float [] outPlaneData = new float[nPlaneRows*nPlaneCols];
        for( int row=0; row<nPlaneRows; row++)
         {
             for( int col=0; col<nPlaneCols; col++)
                 outPlaneData[col+row*nPlaneCols] = applyKernel( inPlaneData, row, col, nPlanes, nPlaneRows,
                                                                 nPlaneCols, rowSize, colSize, kernel);
         }

        return outPlaneData;
    }

    static private float applyKernel( float[][] inPlaneData, int row, int col, int nPlanes, int nPlaneRows,
                                        int nPlaneCols, int rowSize, int colSize, Kernel kernel)
    {
        int u1 = Math.max(col - colSize/2, 0);
        int u2 = Math.min(col + colSize/2, nPlaneCols - 1);
        int v1 = Math.max(row - rowSize/2, 0);
        int v2 = Math.min(row + rowSize/2, nPlaneRows - 1);

        double sum = 0.0;
        for( int plane=0; plane<nPlanes; plane++)
        {
            for( int v=v1; v<=v2; v++)
            {
                for( int u=u1; u<=u2; u++)
                {
                    int kernelRow = u - col + colSize/2; 
                    int kernelCol = v - row + rowSize/2;
                    sum += inPlaneData[plane][u + v*nPlaneCols]*kernel.weights[plane][kernelCol + kernelRow*kernel.width];
                }
            }
        }

        return (float)(sum*kernel.scale);
    }

    static private int getRowSize( int dir, int ySize, int zSize)
    {
        switch(dir)
        {
            case XDIR:
                return zSize;
            case YDIR:
                return zSize;
            case ZDIR:
                return ySize;
            default:
                return 0;
        }
    }

    static private int getColSize( int dir, int xSize, int ySize)
    {
        switch(dir)
        {
            case XDIR:
                return xSize;
            case YDIR:
                return ySize;
            case ZDIR:
                return xSize;
            default:
                return 0;
        }
    }

    /**
     * Compute z = x convolved with y; i.e.,
     *
     *        ifx+lx-1
     *   z[i] = sum  x[j]*y[i-j]  ;  i = ifz,...,ifz+lz-1
     *         j=ifx
     */
    public static void convolve(double[] input, double[] filter, double[] output, int nValues, int firstValue,
                                int lastValue, int windowSize, int windowHalfSize)
    {
        /*
         * Compute output = input convolved with filter; i.e.,
         *
         *          halfWindow + 1
         *   output[i] =   sum    filter[halfWindow + j]*input[i-j],  i = firstValues...1astValue
         *           j=-halfWindow

         * Input:
         * input	     input values
         * firstValue	 sample index of first input value
         * lastValue    sample index of last input value
         * filter       filter values to be convolved with input
         * windowSize	 size of filter array
         * windowCenter center index of filter array
         *
         * Output:
         * ouput		 output values array
         */

        int fMin;
        for(int i = firstValue; i <= lastValue; i++)
        {
            int jMin = i - windowHalfSize;
            int jMax = i + windowHalfSize;

            if(jMin < 0)
            {
                fMin = -jMin;
                jMin = 0;
            }
            else
                fMin = 0;

            if(jMax >= nValues)
                jMax = nValues - 1;

            int f = fMin;
            double sum = 0.0;
            for(int j = jMin; j <= jMax; j++, f++)
                sum += input[j] * filter[f];

            output[i] = sum;
        }
    }
         public static final void convolve(float[] input, double[] filter, float[] output, int nValues, int firstValue,
                                 int lastValue, int windowSize, int windowHalfSize)
     {
         int fMin;
         for(int i = firstValue; i <= lastValue; i++)
         {
             int jMin = i - windowHalfSize;
             int jMax = i + windowHalfSize;

             if(jMin < 0)
             {
                 fMin = -jMin;
                 jMin = 0;
             }
             else
                 fMin = 0;

             if(jMax >= nValues)
                 jMax = nValues - 1;

             int f = fMin;
             double sum = 0.0;
             for(int j = jMin; j <= jMax; j++, f++)
                 sum += input[j] * filter[f];

             output[i] = (float)sum;
         }
     }


    static private double getScale(Kernel kernel)
    {
        double scale = 1;
        double sum = 0.0;
        int planeSize = kernel.width*kernel.height;
        for( int j=0; j<kernel.depth; j++)
        {
            for (int i=0; i<planeSize; i++)
                 sum += kernel.weights[j][i];
        }

        if (sum!=0.0)
             scale = 1.0/sum;

        return scale;
    }

    static public class Kernel
    {
        int width;
        int height;
        int depth;
        double scale;
        float[][] weights;

        public Kernel(int dir, String name, int xSize, int ySize, int zSize)
        {
            setKernelSize( dir, xSize, ySize, zSize);

            if(name.equals(BARTLETT))
            {
                getBartlettWeights();
                scale = getScale(this);
            }
            else if(name.equals(GAUSSIAN))
            {
                getGaussianWeights();
                scale = getScale(this);
            }
            else if(name.equals(LAPLACIAN))
            {
                width  = 3;
                height = 3;
                depth  = 3;
                weights = new float[3][];
                for(int n = 0; n < 3; n++)
                    weights= new float[][] { {-1, 2, -1, 2, -4, 2, -1, 2, -1 },  {2, -4, 2, -4, 8, -4, 2, -4, 2 },  {-1, 2, -1, 2, -4, 2, -1, 2, -1 } };
                scale = 1;
            }
        }

        private void setKernelSize( int dir, int xSize, int ySize, int zSize)
        {
            if( dir==XDIR)
            {
                width  = xSize;
                height = zSize;
                depth  = ySize;
            }
            else if( dir==YDIR)
            {
                width  = ySize;
                height = zSize;
                depth  = xSize;
            }
            else if( dir==ZDIR)
            {
                width  = xSize;
                height = ySize;
                depth  = zSize;
            }
        }

        private void getBartlettWeights()
        {
            weights = new float[depth][width*height];

            float[] vector1 = new float[width];
            float[] vector2 = new float[height];
            float[] vector3 = new float[depth];

            // Initialize first half of Bartlett window
            int i,j,k;
            for( i=0; i<width/2+1; i++)
                vector1[i] = i+1;
            for( j=0; j<height/2+1; j++)
                vector2[j] = j+1;
            for( k=0; k<depth/2+1; k++)
                vector3[k] = k+1;

            // Initialize second half of Bartlett window
            for( i=width/2+1; i<width; i++)
                vector1[i] = width - i;
            for( j=height/2+1; j<height; j++)
                vector2[j] = height - j;
            for( k=depth/2+1; k<depth; k++)
                vector3[k] = depth - k;

            weights = getWeights( vector1, vector2, vector3);
         }

        private void getGaussianWeights()
        {
            weights = new float[depth][width*height];

            float[] vector1 = new float[width];
            float[] vector2 = new float[height];
            float[] vector3 = new float[depth];

            double alpha = 2.0;

            for( int i=-width/2; i<=width/2; i++)
                vector1[i + width/2] = (float)Math.exp(-alpha*((double)i/width)*((double)i/width));
            for( int j=-height/2; j<=height/2; j++)
                vector2[j + height/2] = (float)Math.exp(-alpha*((double)j/height)*((double)j/height));
            for( int k=-depth/2; k<=depth/2; k++)
                vector3[k + depth/2] = (float)Math.exp(-alpha*((double)k/depth)*((double)k/depth));

            weights = getWeights( vector1, vector2, vector3);
        }

        private float[][] getWeights( float[] vector1, float[] vector2, float[] vector3)
        {
            float[][] volumeWeight = new float[depth][width*height];

            if( depth==1)
            {
                volumeWeight[0] = vectorMult(vector1, width, vector2, height);
                return volumeWeight;
            }

            float[] plane1 = vectorMult(vector1, width, vector2, height);
            float[] planeTranspose = transpose2d(plane1, width, height);

            // Compute weights in depth-height planes
            float[] transposeVec = new float[height];
            float[][] volume = new float[width][];
            for( int i=0; i<width; i++)
            {
                for( int j=0; j<height; j++)
                    transposeVec[j] = planeTranspose[j + i*height];
                volume[i] = vectorMult( vector3, depth, transposeVec, height);
            }

            // Sort the weights into height-width planes
            for( int plane=0; plane<depth; plane++)
            {
                for( int row=0; row<height; row++)
                {
                    for( int col=0;  col<width; col++)
                        volumeWeight[plane][col + row*width] = volume[col][plane + row*depth];
                }
            }
            return volumeWeight;
        }

        private float[] vectorMult( float[] vector1, int length1, float[] vector2, int length2)
        {
            float[] productVector = new float[length1*length2];
            for( int i=0; i<length2; i++)
            {
                for( int j=0; j<length1; j++)
                    productVector[j + i*length1] = vector1[j]*vector2[i];
            }
            return productVector;
        }

        private float[] transpose2d( float[] input, int nRow, int nCol)
        {
            float[] transpose = new float[nRow*nCol];
            for( int row = 0; row < nRow; row++)
                for( int col = 0; col < nCol; col++)
                    transpose[row*nCol + col] = input[col*nRow + row];

           return transpose;
        }
     }
 }
