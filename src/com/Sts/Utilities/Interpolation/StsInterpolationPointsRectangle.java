package com.Sts.Utilities.Interpolation;

import com.Sts.Utilities.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: May 28, 2008
 * Time: 7:09:47 PM
 * To change this template use File | Settings | File Templates.
 */
    class StsInterpolationPointsRectangle
    {
        int minX, maxX, minY, maxY;
        int nPoints;
        StsInterpolationPointsRectangle()
        {
            initialize();
        }

        void initialize()
        {
            minX = StsParameters.largeInt;
            maxX = -StsParameters.largeInt;
            minY = StsParameters.largeInt;
            maxY = -StsParameters.largeInt;
            nPoints = 0;
        }

        void addPoint(int x, int y)
        {
            if(nPoints == 0)
            {
                minX = x;
                maxX = x;
                minY = y;
                maxY = y;
            }
            else
            {
                if(x < minX) minX = x;
                else if(x > maxX) maxX = x;
                if(y < minY) minY = y;
                else if(y > maxY) maxY = y;
            }
            nPoints++;
        }

        boolean isInside()
        {
            return minX <= 0 && maxX >= 0 && minY <= 0 && maxY >= 0;
        }

        byte getFillFlag()
        {
            if(isInside())
                return StsParameters.SURF_GAP_FILLED;
            else
                return StsParameters.SURF_BOUNDARY;
        }
    }