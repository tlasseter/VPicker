package com.Sts.UI.Beans;

/**
 * ArrowBeanCustomizer.java  1.00 97/07/09 Merlin Hughes
 *
 * Copyright (c) 1997 Merlin Hughes, All Rights Reserved.
 *
 * Permission to use, copy, modify, and distribute this software
 * for commercial and non-commercial purposes and without fee is
 * hereby granted provided that this copyright notice appears in
 * all copies.
 *
 * http://prominence.com/                         ego@merlin.org
 */

import java.awt.*;
import java.awt.event.*;
import java.beans.*;

public class ArrowBeanCustomizer extends Container implements Customizer, ActionListener {
  ArrowBean left, right, target;

  public ArrowBeanCustomizer () {
    setLayout (new BorderLayout ());;
    try {
      left = (ArrowBean) Beans.instantiate (getClass ().getClassLoader (), "com.Sts.UI.Beans.ArrowBean");
      right = (ArrowBean) Beans.instantiate (getClass ().getClassLoader (), "com.Sts.UI.Beans.ArrowBean");
    } catch (Exception ex) {
      ex.printStackTrace ();
    }
    add ("West", (Component) Beans.getInstanceOf (left, Component.class));
    left.addActionListener (this);
    add ("Center", new Label ("direction", Label.CENTER));
    add ("East", (Component) Beans.getInstanceOf (right, Component.class));
    right.setDirection (ArrowBean.RIGHT);
    right.addActionListener (this);
  }

  public void setObject (Object bean) {
    target = (ArrowBean) bean;
  }

  public void actionPerformed (ActionEvent e) {
    target.setDirection ((e.getSource () == left) ? ArrowBean.LEFT : ArrowBean.RIGHT);
    firePropertyChange ();
  }

  protected PropertyChangeSupport listeners = new PropertyChangeSupport (this);

  public void addPropertyChangeListener (PropertyChangeListener l) {
    listeners.addPropertyChangeListener (l);
  }

  public void removePropertyChangeListener (PropertyChangeListener l) {
    listeners.removePropertyChangeListener (l);
  }

  protected void firePropertyChange () {
    listeners.firePropertyChange ("", null, null);
  }
}
