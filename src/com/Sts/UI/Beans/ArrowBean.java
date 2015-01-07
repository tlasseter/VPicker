package com.Sts.UI.Beans;


import java.awt.*;
import java.awt.event.*;

public class ArrowBean extends Component {
  public static final String LEFT = "left", RIGHT = "right";

  public ArrowBean () {
    enableEvents (AWTEvent.MOUSE_EVENT_MASK);
  }

  protected String direction = LEFT;

  public void setDirection (String d) {
    direction = LEFT.equals (d) ? LEFT : RIGHT;
    repaint ();
  }

  public String getDirection () {
    return direction;
  }
  
  public Dimension getMinimumSize () {
    return getPreferredSize ();
  }

  public Dimension getPreferredSize () {
    return new Dimension (16, 16);
  }

  public Dimension getMaximumSize () {
    return getPreferredSize ();
  }

  protected String od = direction;

  public void update (Graphics g) {
    if (od == direction)
      paint (g);
    else {
      super.update (g);
      od = direction;
    }
  }

  public void paint (Graphics g) {
    int x = getSize ().width, y = getSize ().height;
    if (direction == LEFT) {
      g.setColor ((in && down) ? getBackground ().brighter () :
                  getBackground ().darker ());
      g.drawLine (x - 2, 2, x - 2, y - 2);
      g.drawLine (x - 2, y - 2, 2, (y - 2) / 2);
      g.setColor ((in && down) ? getBackground ().darker () :
                  getBackground().brighter());
      g.drawLine (2, (y - 2) / 2, x - 2, 2);
    } else {
      g.setColor ((in && down) ? getBackground ().brighter () :
                  getBackground ().darker ());
      g.drawLine (2, 2, x - 2, (y - 2) / 2);
      g.drawLine (2, y - 2, x - 2, (y - 2) / 2);
      g.setColor ((in && down) ? getBackground ().darker () :
                  getBackground().brighter());
      g.drawLine (2, 2, 2, y - 2);
    }
  }

  protected static ArrowRoller roller;

  protected synchronized static ArrowRoller createRoller () {
    if (roller == null)
      roller = new ArrowRoller ();
    return roller;
  }

  public synchronized static void destroyRoller () {
    roller.abort ();
    roller = null;
  }

  protected boolean down, in;

  protected void processMouseEvent (MouseEvent e) {
    ArrowRoller roller = this.roller;
    if (roller == null)
      roller = createRoller ();
    switch (e.getID ()) {
      case MouseEvent.MOUSE_PRESSED:
        fireActionEvent ();
        roller.addTarget (this);
        down = in = true;
        repaint ();
        break;
      case MouseEvent.MOUSE_RELEASED:
        roller.removeTarget (this);
        down = in = false;
        repaint ();
        break;
      case MouseEvent.MOUSE_ENTERED:
        in = true;
        if (down) {
          roller.addTarget (this);
          repaint ();
        }
        break;
      case MouseEvent.MOUSE_EXITED:
        in = false;
        if (down) {
          roller.removeTarget (this);
          repaint ();
        }
        break;
    }
    super.processMouseEvent (e);
  }

  protected ActionListener listeners;

  public void addActionListener (ActionListener l) {
    listeners = AWTEventMulticaster.add (l, listeners);
  }

  public void removeActionListener (ActionListener l) {
    listeners = AWTEventMulticaster.remove (l, listeners);
  }

  protected void fireActionEvent () {
    if (listeners != null)
      listeners.actionPerformed (
        new ActionEvent (this, ActionEvent.ACTION_PERFORMED, direction));
  }
}
