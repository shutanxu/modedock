/*
 * molTreeWidget.cpp
 *
 *  Created on: May 31, 2014
 *      Author: stan
 */

#include"molTreeViewer.h"

/*
void
MyTreeWidget::onCustomContextMenuRequested(const QPoint& pos) {
	item = static_cast<MyTreeWidgetItem*>( itemAt(pos) );
	showContextMenu(item, viewport()->mapToGlobal(pos));
	update();
}

void
MyTreeWidget::showContextMenu( QTreeWidgetItem* item, const QPoint& globalPos) {
//	QMenu menu;
//
//	showLineAction = new QAction( tr("Show Line"), this );
//	showLineAction->setCheckable( true );
//	if( item->text(7) == QString::number(1) ){
//		showLineAction->setChecked( true );
//	}else if(  item->text(7) == QString::number(0)  ){
//		showLineAction->setChecked( false );
//	}
//	connect(showLineAction, SIGNAL(toggled(bool)), this, SLOT(inverse_lines()) );
//
//	showStickAction = new QAction( tr("Show Stick"), this );
//	showStickAction->setCheckable(true);
//	if( item->text(8) == QString::number(1) ){
//		showStickAction->setChecked( true );
//	}else if(  item->text(8) == QString::number(0)  ){
//		showStickAction->setChecked( false );
//	}
//	connect(showStickAction, SIGNAL(toggled(bool)), this, SLOT(inverse_sticks()) );
//
//	showVdwBallAction = new QAction( tr("Show VDW Ball"), this );
//	showVdwBallAction->setCheckable(true);
//	if( item->text(9) == QString::number(1) ){
//		showVdwBallAction->setChecked( true );
//	}else if(  item->text(9) == QString::number(0)  ){
//		showVdwBallAction->setChecked( false );
//	}
//	connect(showVdwBallAction, SIGNAL(toggled(bool)), this, SLOT(inverse_vdwBall()) );
//
//	showStickBallAction = new QAction( tr("Show Stick Ball"), this );
//	showStickBallAction->setCheckable(true);
//	if( item->text(10) == QString::number(1) ){
//		showStickBallAction->setChecked( true );
//	}else if(  item->text(10) == QString::number(0)  ){
//		showStickBallAction->setChecked( false );
//	}
//	connect(showStickBallAction, SIGNAL(toggled(bool)), this, SLOT(inverse_stickBall()) );
//
//	menu.addAction( showLineAction );
//	menu.addAction( showStickAction );
//	menu.addAction( showVdwBallAction );
//	menu.addAction( showStickBallAction );
//	menu.exec(globalPos);
}

void
MyTreeWidget::inverse_lines(){
	if( item->text(7) == QString::number(1) ){
		item->setText(7, QString::number(0) );

		//update children
		for( int i=0; i<item->childCount(); i++ ){
			item->child( i )->setText(7, QString::number(0) );
			for( int j=0; j<item->child(i)->childCount(); j++ ){
				item->child( i )->child(j)->setText(7, QString::number(0) );
				for( int k=0; k<item->child(i)->child(j)->childCount(); k++ ){
					item->child( i )->child(j)->child(k)->setText(7, QString::number(0) );
				}
			}
		}

	}else if(  item->text(7) == QString::number(0)  ){
		item->setText(7, QString::number(1) );

		//update children
		for( int i=0; i<item->childCount(); i++ ){
			item->child( i )->setText(7, QString::number(1) );
			for( int j=0; j<item->child(i)->childCount(); j++ ){
				item->child( i )->child(j)->setText(7, QString::number(1) );
				for( int k=0; k<item->child(i)->child(j)->childCount(); k++ ){
					item->child( i )->child(j)->child(k)->setText(7, QString::number(1) );
				}
			}
		}
	}
	update();
}

void
MyTreeWidget::inverse_sticks(){
	if( item->text(8) == QString::number(1) ){
		item->setText(8, QString::number(0) );
		//update children
		for( int i=0; i<item->childCount(); i++ ){
			item->child( i )->setText(8, QString::number(0) );
			for( int j=0; j<item->child(i)->childCount(); j++ ){
				item->child( i )->child(j)->setText(8, QString::number(0) );
				for( int k=0; k<item->child(i)->child(j)->childCount(); k++ ){
					item->child( i )->child(j)->child(k)->setText(8, QString::number(0) );
				}
			}
		}
	}else if(  item->text(8) == QString::number(0)  ){
		item->setText(8, QString::number(1) );
		//update children
		for( int i=0; i<item->childCount(); i++ ){
			item->child( i )->setText(8, QString::number(1) );
			for( int j=0; j<item->child(i)->childCount(); j++ ){
				item->child( i )->child(j)->setText(8, QString::number(1) );
				for( int k=0; k<item->child(i)->child(j)->childCount(); k++ ){
					item->child( i )->child(j)->child(k)->setText(8, QString::number(1) );
				}
			}
		}
	}
	update();
}

void
MyTreeWidget::inverse_vdwBall(){
	if( item->text(9) == QString::number(1) ){
		item->setText(9, QString::number(0) );
		//update children
		for( int i=0; i<item->childCount(); i++ ){
			item->child( i )->setText(9, QString::number(0) );
			for( int j=0; j<item->child(i)->childCount(); j++ ){
				item->child( i )->child(j)->setText(9, QString::number(0) );
				for( int k=0; k<item->child(i)->child(j)->childCount(); k++ ){
					item->child( i )->child(j)->child(k)->setText(9, QString::number(0) );
				}
			}
		}
	}else if(  item->text(9) == QString::number(0)  ){
		item->setText(9, QString::number(1) );
		//update children
		for( int i=0; i<item->childCount(); i++ ){
			item->child( i )->setText(9, QString::number(1) );
			for( int j=0; j<item->child(i)->childCount(); j++ ){
				item->child( i )->child(j)->setText(9, QString::number(1) );
				for( int k=0; k<item->child(i)->child(j)->childCount(); k++ ){
					item->child( i )->child(j)->child(k)->setText(9, QString::number(1) );
				}
			}
		}
	}
	update();
}

void
MyTreeWidget::inverse_stickBall(){
	if( item->text(10) == QString::number(1) ){
		item->setText(10, QString::number(0) );
		//update children
		for( int i=0; i<item->childCount(); i++ ){
			item->child( i )->setText(10, QString::number(0) );
			for( int j=0; j<item->child(i)->childCount(); j++ ){
				item->child( i )->child(j)->setText(10, QString::number(0) );
				for( int k=0; k<item->child(i)->child(j)->childCount(); k++ ){
					item->child( i )->child(j)->child(k)->setText(10, QString::number(0) );
				}
			}
		}
	}else if(  item->text(10) == QString::number(0)  ){
		item->setText(10, QString::number(1) );
		//update children
		for( int i=0; i<item->childCount(); i++ ){
			item->child( i )->setText(10, QString::number(1) );
			for( int j=0; j<item->child(i)->childCount(); j++ ){
				item->child( i )->child(j)->setText(10, QString::number(1) );
				for( int k=0; k<item->child(i)->child(j)->childCount(); k++ ){
					item->child( i )->child(j)->child(k)->setText(10, QString::number(1) );
				}
			}
		}
	}
	update();
}

/////////////////////////////////////////////////////////////////////////////////

MolTreeViewer::MolTreeViewer( const vector<Molecular>& allMol, QWidget *parent ){

    for( size_t i=0; i<allMol.size(); i++ ){
    	molVec.push_back( allMol[i] );
    }

	treeWidget = new MyTreeWidget;
	treeWidget->setColumnCount( 1 );
	treeWidget->setHeaderLabels( QStringList() << tr("Molecular") );
//    treeWidget->header()->setResizeMode(0, QHeaderView::Stretch);

    addMol();

//    QVBoxLayout *mainLayout = new QVBoxLayout;
//    mainLayout->addWidget(treeWidget);
//    mainLayout->setMargin(0);
//    setLayout(mainLayout);
//	setFocus();
}

void
MolTreeViewer::setMolVec( const vector<Molecular>& mVec ){
	molVec.clear();
	molVec = mVec;
    addMol();
}

void
MolTreeViewer::addMol ( ){

	allWidgetItemVec.clear();
	treeWidget->clear();

	for( size_t i=0; i<molVec.size(); i++ ){
		Molecular						mol = molVec[i];
		MyTreeWidgetItem		*itm = new MyTreeWidgetItem( treeWidget );
		string							molName = molVec[i].get_name();
		itm->setCheckState(0, Qt::Checked);
		itm->setText( 0, QString( molName.c_str() ) );		// molecular Index;
		itm->setText( 1, QString::number(i) );						// molecular Index;
		itm->setText( 2, QString("NULL") );							// chain name						 indicate the chain father
		itm->setText( 3, QString::number(-1) );					// subChain type;
		itm->setText( 4, QString::number(-1) );					// residue index
		itm->setText( 5, QString::number(-1) );					// atom index
		itm->setText( 6, QString("NULL") );							// atom name
		itm->setText( 7, QString::number( 1 ) );					// indicate whether to show line model, 1 to show, 0 not to show
		itm->setText( 8, QString::number( 0 ) );					// indicate whether to show stick model, 1 to show, 0 not to show
		itm->setText( 9, QString::number( 0 ) );					// indicate whether to show VDW ball model, 1 to show, 0 not to show
		itm->setText( 10, QString::number( 0 ) );				// indicate whether to show stick ball model, 1 to show, 0 not to show
		allWidgetItemVec.push_back( itm );
		for( size_t j=0; j<mol.get_chainVec().size(); j++ ){
			addChain( itm, i, mol.get_chainVec()[j] );
		}
	}
}

void
MolTreeViewer::addChain( MyTreeWidgetItem *parent,
													const size_t& moleIndex,const Chain& ch ){
	if( ch.get_subChainVec().size() == 0 ){
		return;
	}
	MyTreeWidgetItem *itm = new MyTreeWidgetItem();
	itm->setCheckState(0, Qt::Checked);
	if( !ch.get_name().empty() ){
		string					chName = ch.get_name();

		itm->setText( 0, QString(ch.get_name().c_str()) );					//
		itm->setText( 1, QString::number( moleIndex ) );						// molecular Index;
		itm->setText( 2, QString(  chName.c_str()  ) );							// chain name						 indicate the chain father
		itm->setText( 3, QString::number(-1) );										// subChain type;
		itm->setText( 4, QString::number(-1) );										// residue index
		itm->setText( 5, QString::number(-1) );										// atom index
		itm->setText( 6, QString("NULL") );												// atom name
		itm->setText( 7, QString::number( 1 ) );										// indicate whether to show line model, 1 to show, 0 not to show
		itm->setText( 8, QString::number( 0 ) );										// indicate whether to show stick model, 1 to show, 0 not to show
		itm->setText( 9, QString::number( 0 ) );										// indicate whether to show VDW ball model, 1 to show, 0 not to show
		itm->setText( 10, QString::number( 0 ) );									// indicate whether to show stick ball model, 1 to show, 0 not to show
	}
	allWidgetItemVec.push_back( itm );

	for( size_t i=0; i<ch.get_subChainVec().size(); i++ ){
		addSubChain( itm, moleIndex, ch, ch.get_subChainVec()[i] );
	}
	parent->addChild( itm );
}

void
MolTreeViewer::addSubChain( MyTreeWidgetItem *parent,
		const size_t& moleIndex,
		const Chain& ch,
		const SubChain& sub ){
	if( sub.get_ResidueVec().size() == 0 ){
		return;
	}
	string chName = sub.get_ResidueVec().front().get_chainName();
	MyTreeWidgetItem *itm = new MyTreeWidgetItem();
	itm->setCheckState(0, Qt::Checked);

	if( !sub.get_ResidueVec().empty() ){
		int typeIndex = sub.get_molType();
		string		type;
		if( typeIndex == 0 ){
			type = string("PROTEIN");
		}else if( typeIndex == 1 ){
			type = string("DNA");
		}else if( typeIndex == 2 ){
			type = string("RNA");
		}else if( typeIndex == 3 ){
			type = string("LIGAND");
		}else if( typeIndex == 4 ){
			type = string("HOH");
		}
		itm->setText( 0,  QString(  type.c_str()  ) );
		itm->setText( 1, QString::number( moleIndex ) );	// molecular Index;
		itm->setText( 2, QString(  chName.c_str()  ) );		// chain name						 indicate the chain father
		itm->setText( 3, QString::number( sub.get_molType() ) );							// subChain type;
		itm->setText( 4, QString::number(-1) );					// residue index
		itm->setText( 5, QString::number(-1) );					// atom index
		itm->setText( 6, QString("NULL") );							// atom name
		itm->setText( 7, QString::number( 1 ) );					// indicate whether to show line model, 1 to show, 0 not to show
		itm->setText( 8, QString::number( 0 ) );					// indicate whether to show stick model, 1 to show, 0 not to show
		itm->setText( 9, QString::number( 0 ) );					// indicate whether to show VDW ball model, 1 to show, 0 not to show
		itm->setText( 10, QString::number( 0 ) );				// indicate whether to show stick ball model, 1 to show, 0 not to show
	}

	allWidgetItemVec.push_back( itm );
	vector<Residue>			resVec = sub.get_ResidueVec();
	for( size_t i=0; i<resVec.size(); i++ ){
		addResidue( itm, moleIndex, ch, sub, resVec[i] );
	}
	parent->addChild( itm );
}

void
MolTreeViewer::addResidue( MyTreeWidgetItem *parent,
		const size_t& moleIndex,
		const Chain& ch,
		const SubChain& sub ,
		const Residue& res){
	string						chName = res.get_chainName();
	vector<Atom>		atVec = res.get_atomVec();
	MyTreeWidgetItem *itm = new MyTreeWidgetItem();
	itm->setCheckState(0, Qt::Checked);
	stringstream			ss;
	ss<<int( res.get_index() );
	string						residName = ss.str() +" "+res.get_name();

	itm->setText( 0,  QString::fromStdString( residName ) );
	itm->setText( 1, QString::number( moleIndex ) );	// molecular Index;
	itm->setText( 2, QString(  chName.c_str()  ) );		// chain name						 indicate the chain father
	itm->setText( 3, QString::number( sub.get_molType() ) );							// subChain type;
	itm->setText( 4, QString::number( res.get_index()) );					// residue index
	itm->setText( 5, QString::number(-1) );					// atom index
	itm->setText( 6, QString("NULL") );							// atom name
	itm->setText( 7, QString::number( 1 ) );					// indicate whether to show line model, 1 to show, 0 not to show
	itm->setText( 8, QString::number( 0 ) );					// indicate whether to show stick model, 1 to show, 0 not to show
	itm->setText( 9, QString::number( 0 ) );					// indicate whether to show VDW ball model, 1 to show, 0 not to show
	itm->setText( 10, QString::number( 0 ) );				// indicate whether to show stick ball model, 1 to show, 0 not to show

	allWidgetItemVec.push_back( itm );
	for( size_t i=0; i<atVec.size(); i++ ){
		addAtom( itm, moleIndex, ch, sub, res, atVec[i] );
	}
	parent->addChild( itm );
}

void
MolTreeViewer::addAtom(  MyTreeWidgetItem *parent,
		const size_t& moleIndex,
		const Chain& ch,
		const SubChain& sub ,
		const Residue& res,
		const Atom& at  ){
	MyTreeWidgetItem *itm = new MyTreeWidgetItem();
	itm->setCheckState(0, Qt::Checked);
	stringstream		s1;
	stringstream		ss;
	string					s;
	s1.clear();
	s1.width(6);
	s1<<at.get_index();

	s1.clear();
	s1.width(5);
	s1<<at.get_name();

	s1.clear();
	s1.width(16);
	s1<<at.get_coord().x ;

	s1.clear();
	s1.width(16);
	s1<<at.get_coord().y;

	s1.clear();
	s1.width(16);
	s1<<at.get_coord().z;
	s =s1.str();

	itm->set_atom( at );

	itm->setText( 0,  QString::fromStdString( s ) );
	itm->setText( 1, QString::number( moleIndex ) );	// molecular Index;
	itm->setText( 2, QString(  ch.get_name().c_str()  ) );		// chain name						 indicate the chain father
	itm->setText( 3, QString::number( sub.get_molType() ) );							// subChain type;
	itm->setText( 4, QString::number( res.get_index()) );					// residue index
	itm->setText( 5, QString::number(  at.get_index() ) );					// atom index
	itm->setText( 6, QString::fromStdString(at.get_name()) );		// atom name
	itm->setText( 7, QString::number( 1 ) );					// indicate whether to show line model, 1 to show, 0 not to show
	itm->setText( 8, QString::number( 0 ) );					// indicate whether to show stick model, 1 to show, 0 not to show
	itm->setText( 9, QString::number( 0 ) );					// indicate whether to show VDW ball model, 1 to show, 0 not to show
	itm->setText( 10, QString::number( 0 ) );				// indicate whether to show stick ball model, 1 to show, 0 not to show

	allWidgetItemVec.push_back( itm );

	parent->addChild( itm );
}

void
MolTreeViewer::updateChildItemState( MyTreeWidgetItem* item ){

	// molecular, if Uncheck molecular checkbox, then all chains under it are Unchecked
	if( item->checkState(0) == Qt::Unchecked  && item->text(2) == QString("NULL")  ){
		for( size_t i=0; i<allWidgetItemVec.size(); i++ ){
			if( allWidgetItemVec[i]->text(1) == item->text(1) ){
				allWidgetItemVec[i]->setData(0, Qt::CheckStateRole, Qt::Unchecked);
			}
		}
	}else	if( item->checkState(0) == Qt::Checked  && item->text(2) == QString("NULL") ){
		for( size_t i=0; i<allWidgetItemVec.size(); i++ ){
			if( allWidgetItemVec[i]->text(1) == item->text(1) ){
				allWidgetItemVec[i]->setData(0, Qt::CheckStateRole, Qt::Checked);
			}
		}
	}

	// chain, if Uncheck chain checkbox, then all sub chain under it are Unchecked
	if( item->checkState(0) == Qt::Unchecked  &&
			item->text(2) != QString("NULL")  &&
			item->text(3) == QString::number(-1)  ){
		for( size_t i=0; i<allWidgetItemVec.size(); i++ ){
			if( allWidgetItemVec[i]->text(2) == item->text(2) &&
					allWidgetItemVec[i]->text(1) == item->text(1)  ){
				allWidgetItemVec[i]->setData(0, Qt::CheckStateRole, Qt::Unchecked);
			}
		}
	}else	if( item->checkState(0) == Qt::Checked &&
			item->text(2) != QString("NULL")  &&
			item->text(3) == QString::number(-1) ){
		for( size_t i=0; i<allWidgetItemVec.size(); i++ ){
			if( allWidgetItemVec[i]->text(2) == item->text(2) &&
					allWidgetItemVec[i]->text(1) == item->text(1)  ){
				allWidgetItemVec[i]->setData(0, Qt::CheckStateRole, Qt::Checked);
			}
		}
	}

	// chain, if Uncheck sub chain checkbox, then all residues under it are Unchecked
	if( item->checkState(0) == Qt::Unchecked  &&
			item->text(2) != QString("NULL")  &&
			item->text(3) != QString::number(-1) &&
			item->text(4) == QString::number(-1)  ){
		for( size_t i=0; i<allWidgetItemVec.size(); i++ ){
			if( allWidgetItemVec[i]->text(2) == item->text(2) &&
					allWidgetItemVec[i]->text(1) == item->text(1) &&
					allWidgetItemVec[i]->text(3) == item->text(3)  ){
				allWidgetItemVec[i]->setData(0, Qt::CheckStateRole, Qt::Unchecked);
			}
		}
	}else	if( item->checkState(0) == Qt::Checked &&
			item->text(2) != QString("NULL")  &&
			item->text(3) != QString::number(-1)  &&
			item->text(4) == QString::number(-1) ){
		for( size_t i=0; i<allWidgetItemVec.size(); i++ ){
			if( allWidgetItemVec[i]->text(2) == item->text(2) &&
					allWidgetItemVec[i]->text(1) == item->text(1) &&
					allWidgetItemVec[i]->text(3) == item->text(3)  ){
				allWidgetItemVec[i]->setData(0, Qt::CheckStateRole, Qt::Checked);
			}
		}
	}

	// chain, if Uncheck residues checkbox, then all atoms under it are Unchecked
	if( item->checkState(0) == Qt::Unchecked  &&
			item->text(2) != QString("NULL")  &&
			item->text(3) != QString::number(-1) &&
			item->text(4) != QString::number(-1)  &&
			item->text(5) == QString::number(-1)  ){
		for( size_t i=0; i<allWidgetItemVec.size(); i++ ){
			if(allWidgetItemVec[i]->text(4) == item->text(4) &&
					allWidgetItemVec[i]->text(3) == item->text(3) &&
					allWidgetItemVec[i]->text(2) == item->text(2) &&
					allWidgetItemVec[i]->text(1) == item->text(1)   ){
				allWidgetItemVec[i]->setData(0, Qt::CheckStateRole, Qt::Unchecked);
			}
		}
	}else	if( item->checkState(0) == Qt::Checked &&
			item->text(2) != QString("NULL")  &&
			item->text(3) != QString::number(-1)  &&
			item->text(4) != QString::number(-1)  &&
			item->text(5) == QString::number(-1) ){
		for( size_t i=0; i<allWidgetItemVec.size(); i++ ){
			if( allWidgetItemVec[i]->text(4) == item->text(4) &&
					allWidgetItemVec[i]->text(3) == item->text(3) &&
					allWidgetItemVec[i]->text(2) == item->text(2) &&
					allWidgetItemVec[i]->text(1) == item->text(1)   ){
				allWidgetItemVec[i]->setData(0, Qt::CheckStateRole, Qt::Checked);
			}
		}
	}
}

vector< vector< pair<Atom, Atom> > >
MolTreeViewer::getVisibleLineModeAtomPairVec(){
	vector< vector<pair<Atom, Atom> > >			atomPairVec;
	for( size_t j=0; j<molVec.size(); j++ ){
		clock_t			t1=clock();
		vector< Atom >												visibleAtomVec;			// pair<reisdIndex, atomName>
		vector<Atom>												molAtomVec = molVec[j].get_atomVec() ;

		size_t																widgeSize = allWidgetItemVec.size();
		for( size_t i=0; i<widgeSize; i++ ){
			if(     allWidgetItemVec[i]->checkState(0) == Qt::Checked &&
					allWidgetItemVec[i]->text(1) ==QString::number( j )&&
					allWidgetItemVec[i]->text(2) != QString( "NULL" ) &&
					allWidgetItemVec[i]->text(3) != QString::number(-1) &&
					allWidgetItemVec[i]->text(4) != QString::number( -1 ) &&
					allWidgetItemVec[i]->text(5) != QString::number( -1 ) &&
					allWidgetItemVec[i]->text(7) == QString::number( 1 )
			){
				Atom										atom = allWidgetItemVec[i]->get_atom();
				visibleAtomVec.push_back( atom );
			}
		}

		vector< pair<Atom, Atom> >	atomPairOneMol = getBondAtomPairs( molVec[j], visibleAtomVec );

		atomPairVec.push_back( atomPairOneMol );
	}

	return 										atomPairVec;
}

vector< vector< pair<Atom, Atom> > >
MolTreeViewer::getVisibleStickModeAtomPairVec(){
	vector< vector<pair<Atom, Atom> > >			atomPairVec;
	for( size_t j=0; j<molVec.size(); j++ ){
		vector< Atom >												visibleAtomVec;			// pair<reisdIndex, atomName>
		vector<Atom>												molAtomVec = molVec[j].get_atomVec() ;
		for( size_t i=0; i<allWidgetItemVec.size(); i++ ){
			if(    allWidgetItemVec[i]->checkState(0) == Qt::Checked &&
					allWidgetItemVec[i]->text(1) ==QString::number( j )&&
					allWidgetItemVec[i]->text(2) != QString("NULL") &&
					allWidgetItemVec[i]->text(3) != QString::number( -1 ) &&
					allWidgetItemVec[i]->text(4) != QString::number( -1 ) &&
					allWidgetItemVec[i]->text(5) != QString::number( -1 ) &&
					allWidgetItemVec[i]->text(8) == QString::number( 1 ) ){
//				Atom										atom = searchAtom( molAtomVec,
//						string2int( allWidgetItemVec[i]->text(5).toStdString() ),
//						allWidgetItemVec[i]->text(6).toStdString() ,
//						string2int(allWidgetItemVec[i]->text(4).toStdString() ),
//						allWidgetItemVec[i]->text(2).toStdString() 	);
				Atom										atom = allWidgetItemVec[i]->get_atom();

				visibleAtomVec.push_back( atom );
			}
		}
		vector< pair<Atom, Atom> >	atomPairOneMol = getBondAtomPairs( molVec[j], visibleAtomVec );
		atomPairVec.push_back( atomPairOneMol );
	}
	return 										atomPairVec;
}

vector< vector< Atom > >
MolTreeViewer::getVisibleVDWballModeAtomPairVec(){
	vector< vector<Atom > >			atomPairVec;
	clock_t			t1 = clock() ;
	for( size_t j=0; j<molVec.size(); j++ ){
		vector<size_t>						visibleAtomVec;
		vector<string>						visibleAtomNameVec;
		vector<size_t>						visibleAtomResidIndexVec;
		for( size_t i=0; i<allWidgetItemVec.size(); i++ ){
			if(   allWidgetItemVec[i]->checkState(0) == Qt::Checked &&
					allWidgetItemVec[i]->text(1) ==QString::number( j )&&
					allWidgetItemVec[i]->text(2) != QString("NULL") &&
					allWidgetItemVec[i]->text(3) != QString::number( -1 ) &&
					allWidgetItemVec[i]->text(4) != QString::number( -1 ) &&
					allWidgetItemVec[i]->text(5) != QString::number( -1 ) &&
					allWidgetItemVec[i]->text(9) == QString::number( 1 ) ){

				visibleAtomVec.push_back( string2int( allWidgetItemVec[i]->text(5).toStdString() ) );
				visibleAtomNameVec.push_back( allWidgetItemVec[i]->text(6).toStdString()  );
				visibleAtomResidIndexVec.push_back( string2int( allWidgetItemVec[i]->text(4).toStdString() ) );
			}
		}

		vector< Atom >	atomPairOneMol = getSubAtomVec( molVec[j], visibleAtomNameVec,
				visibleAtomVec,  visibleAtomResidIndexVec );
		atomPairVec.push_back( atomPairOneMol );
	}
	return 										atomPairVec;
}

vector< vector< pair<Atom, Atom> > >
MolTreeViewer::getVisibleStickBallModeAtomPairVec(){
	vector< vector<pair<Atom, Atom> > >			atomPairVec;
//	vector< vector<> >				allVisibleAtomVec;
	for( size_t j=0; j<molVec.size(); j++ ){
		vector< pair<size_t, string> >						visibleAtomVec;			// pair<reisdIndex, atomName>
		for( size_t i=0; i<allWidgetItemVec.size(); i++ ){
			if(   allWidgetItemVec[i]->checkState(0) == Qt::Checked &&
					allWidgetItemVec[i]->text(1) ==QString::number( j )&&
					allWidgetItemVec[i]->text(2) != QString("NULL") &&
					allWidgetItemVec[i]->text(3) != QString::number( -1 ) &&
					allWidgetItemVec[i]->text(4) != QString::number( -1 ) &&
					allWidgetItemVec[i]->text(5) != QString::number( -1 ) &&
					allWidgetItemVec[i]->text(10) == QString::number( 1 )  ){
				pair<size_t, string >				atom= make_pair<size_t, string >(string2int( allWidgetItemVec[i]->text(4).toStdString() ),
																													allWidgetItemVec[i]->text(6).toStdString()  );
//				Atom										atom = allWidgetItemVec[i]->get_atom();

				visibleAtomVec.push_back( atom );

			}
		}
		vector< pair<Atom, Atom> >	atomPairOneMol = getBondAtomPairs( molVec[j], visibleAtomVec );
		atomPairVec.push_back( atomPairOneMol );
	}
	return 										atomPairVec;
}
*/

//----------------------------------------------------------------------------
MolTreeViewer::MolTreeViewer(){
    createActions();
}

void
MolTreeViewer::createActions(){
//    for( size_t i=0; i<molRow.size(); i++ ){
//        connect( molRow[0], SIGNAL(itemChanged()), this, SLOT( update() ) );
//    }
//    connect( item, SIGNAL(triggered()), this, SLOT( update() ) );
}

void
MolTreeViewer::molChanged(QStandardItem *it){

    if( it->checkState() == Qt::Unchecked ){
        int index = it->row();
        visibleMolIndex[index] = 0;
    }else{
        int index = it->row();
        visibleMolIndex[index] = 1;
    }
    emitMolChanged( );
}

vector<Molecular>
MolTreeViewer::getVisibleMols(){
    vector<Molecular> mols ;
    for( size_t i=0; i<visibleMolIndex.size(); i++ ){
        if(visibleMolIndex[i] != 0){
            mols.push_back( molVec[i] );
        }
    }
    return mols;
}

void
MolTreeViewer::setMolVec( const vector<Molecular> mol ){
	molVec=mol;
    visibleMolIndex.clear();
    standardModel = new QStandardItemModel;
	QStandardItem *rootNode=standardModel->invisibleRootItem();

	for( size_t i=0; i<molVec.size(); i++ ){
        visibleMolIndex.push_back(1);
		QStandardItem *molItem = new QStandardItem( molVec[i].get_name().c_str() );
        molItem->setCheckable(true);
        molItem->setCheckState(Qt::Checked);
        vector<Residue> resVec=molVec[i].get_residueVec();
        for( size_t j=0; j<resVec.size(); j++ ){
            ostringstream resI;
            resI<<resVec[j].get_index();

            string res = resI.str() + string("  ") + resVec[j].get_name();
            QStandardItem *resItem = new QStandardItem( res.c_str());

            vector<Atom> atVec=resVec[j].get_atomVec();
            for( size_t k=0; k<atVec.size(); k++ ){
                QList< QStandardItem* > atItem;

                ostringstream convert;
                convert<<atVec[k].get_index();

                string s=atVec[k].get_name()+convert.str();
                QStandardItem *item=new QStandardItem( s.c_str() );
                item->setEditable(false);
                atItem <<item;
                resItem->appendRow( atItem );
            }
            molItem->appendRow( resItem );
        }
		rootNode->appendRow(molItem);
	}
    connect(standardModel, SIGNAL(itemChanged(QStandardItem*)), this, SLOT(molChanged(QStandardItem*)) ) ;
	setModel( standardModel );
}
