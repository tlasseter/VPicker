package com.Sts.Actions;

import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public class StsMallMap extends JDialog
{
    StsJPanel panel = StsJPanel.addInsets();
    JPanel map = new JPanel();
    StsBooleanFieldBean userChkbox = new StsBooleanFieldBean(this, "users", "Users");
    StsBooleanFieldBean northArrowChkbox = new StsBooleanFieldBean(this, "northArrow", "North Arrow");
    StsBooleanFieldBean projectBoundsChkbox = new StsBooleanFieldBean(this, "projectBounds", "Project Bounds");
    StsBooleanFieldBean objectBoundsChkbox = new StsBooleanFieldBean(this, "objectBounds", "Object Bounds");
    StsBooleanFieldBean fovChkbox = new StsBooleanFieldBean(this, "fieldOfView", "Field of View");

    StsModel model = null;
    Graphics mapGraphics = null;
    boolean users = true;
	boolean displayNorthArrow = true;
	boolean displayProjectBounds = true;
	boolean displayObjectBounds = true;
	boolean displayFieldOfView = true;

    public StsMallMap(Frame frame, String title, boolean modal)
	{
        super(frame, title, modal);
        try
		{
            jbInit();
            pack();
        }
        catch(Exception ex)
		{
            ex.printStackTrace();
        }
    }

    public StsMallMap()
	{
        this(null, "", false);
    }

    public void init(StsModel model)
	{
        this.model = model;
        mapGraphics = map.getGraphics();
        initDisplay();
        this.setVisible(true);
        return;
    }

    public void initDisplay()
	{
        this.repaint();
        return;
    }

    public void paint(Graphics g)
    {
        float[] point = new float[2];
        Color color = new Color(0,255,100);
        super.paint(g);

        mapGraphics.setColor(color);
        mapGraphics.drawRect(this.map.getX()+2, this.map.getY()+2, this.map.getWidth()-4, this.map.getHeight()-4);
        mapGraphics.drawString("Map Here", this.map.getX() + this.map.getWidth()/2 - 20, this.map.getHeight()/2);

        // Users

        // North Arrow

        // Project Boundary
		if(displayProjectBounds)
		{
            StsProject project = model.getProject();
            point[0] = project.getXMin();
			point[1] = project.getXMax();
			int[] p1 = convertPoint(point);
			point[1] = project.getYMin();
			point[0] = project.getYMax();
			int[] p2 = convertPoint(point);
		}
        // Object Boundaries

        // Field of View
    }

    private int[] convertPoint(float[] point)
	{
        int[] ipoint = new int [2];
        return ipoint;
    }

    private void jbInit() throws Exception
	{
        map.setBackground(Color.black);
		map.setForeground(Color.gray);
        map.setBorder(BorderFactory.createLineBorder(Color.black));
		map.setMinimumSize(new Dimension(150, 150));
        getContentPane().add(panel);
        panel.add(map);
		panel.gbc.anchor = GridBagConstraints.WEST;
        panel.add(fovChkbox);
        panel.add(userChkbox);
        panel.add(northArrowChkbox);
        panel.add(projectBoundsChkbox);
        panel.add(objectBoundsChkbox);
    }

	public void setUsers(boolean users) { this.users = users; }
	public boolean getUsers() { return users; }
	public void setProjectBounds(boolean bounds) { displayProjectBounds = bounds; }
	public boolean getProjectBounds() { return displayProjectBounds; }
	public void setObjectBounds(boolean bounds) { displayObjectBounds = bounds; }
	public boolean getObjectBounds() { return displayObjectBounds; }
	public void setFieldOfView(boolean fov) { displayFieldOfView = fov; }
	public boolean getFieldOfView() { return displayFieldOfView; }

	public static void main(String[] args)
	{
		StsMallMap mallMap = new StsMallMap(null, "Mall Map Test", true);
		StsModel model = new StsModel();
		mallMap.init(model);
		mallMap.setVisible(true);
	}
}
