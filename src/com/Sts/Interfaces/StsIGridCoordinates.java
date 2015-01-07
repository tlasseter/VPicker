
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Stuart Jackson
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Interfaces;


public interface StsIGridCoordinates
{
                                          // Inline
    float getRowLabelInc();               // Line Increment
    int getNRows();                       // Number of Lines
    int getRowMin();                      // Min Line Index
    int getRowMax();                      // Max Line Index
    float getRowCoor(float y);            // Get Line Index from Y Coordinate

    double getYOrigin();                  // True Y Origin
    float getYMin();                      // Relative Y Minimum
    float getYMax();                      // Relative Y Maximum
    float getYInc();                      // Distance Between Lines
    float getYSize();                     // Distance Between Min and Max Lines
    float getYCoor(float rowF);           // Get Relative Y From Line Index

                                          // Crossline
    float getColLabelInc();               // Crossline Increment
    int getNCols();                       // Number of Crosslines
    int getColMin();                      // Min Crossline Index
    int getColMax();                      // Max Crossline Index
    float getColCoor(float x);            // Get Crossline Index from X Coordinate

    double getXOrigin();                  // True X Origin
    float getXMin();                      // Relative X Minimum
    float getXMax();                      // Relative X Maximum
    float getXInc();                      // Distance Between Crosslines
    float getXSize();                     // Distance Between Min and Max Crosslines
    float getXCoor(float colF);           // Get Relative X From Crossline Index

                                          // Slice
    int getNSlices();                     // Number of Slices
    int getSliceMin();                    // Min Slice Index
    int getSliceMax();                    // Max Slice Index
    float getSliceCoor(float z);          // Get Slice Index from Z Coordinate

    float getZMin();                      // True Z Min
    float getZMax();                      // True Z Max
    float getZInc();                      // Z (time/depth) Distance Between Slices
    float getZSize();                     // Z (time/depth) Distance Between Min and Max Slices
    float getZCoor(float z);              // Get Z From Slice Index

    float getAngle();                     // Data Rotation Angle

    float getLabelFromCoor(int dir, float coor);    // Get Label from Relative Coordinate
    float getIndexFromLabel(int dir, float label);  // Get Index from Label
    float getLabelFromIndex(int dir, float index);  // Get Label from Index
}
