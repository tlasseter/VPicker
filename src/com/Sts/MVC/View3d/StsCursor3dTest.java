package com.Sts.MVC.View3d;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jul 3, 2008
 * Time: 11:58:12 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsCursor3dTest
{
    public float xMin = -50.0f;
    public float xMax = 50.0f;
    public float xInc = 1.0f;
    public float yMin = -100.0f;
    public float yMax = 100.0f;
    public float yInc = 1.0f;
    public float zMin = -200.0f;
    public float zMax = 200.0f;
    public float zInc = 1.0f;

    public StsCursor3dTest()
    {
    }

    public float getCurrentDirCoordinate(int dir)
    {
        if (dir == StsCursor3d.NONE) return 0.0f;
        else return 50.0f;
    }
}
