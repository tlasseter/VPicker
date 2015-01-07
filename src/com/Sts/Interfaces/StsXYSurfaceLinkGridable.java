package com.Sts.Interfaces;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Aug 12, 2009
 * Time: 4:17:09 PM
 * To change this template use File | Settings | File Templates.
 */

/** Objects implementing this interface provide both an XY orthogonal grid,
 *  but also links which define connections between points on a grid.
 *  Currently used in the construction of a TriangleStrip.
 *  If a link between two points in adjacent rows or columns, doesn't exist,
 *  then this edge on a TStrip doesn't exist.
 */
public interface StsXYSurfaceLinkGridable extends StsXYSurfaceGridable
{
    public boolean hasRowLink(int row, int col);
    public boolean hasColLink(int row, int col);
}
