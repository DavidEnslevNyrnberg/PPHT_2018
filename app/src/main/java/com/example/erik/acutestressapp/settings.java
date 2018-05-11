package com.example.erik.acutestressapp;

import android.os.Bundle;
import android.support.v7.app.AppCompatActivity;
import android.view.View;
import android.widget.CheckBox;

public class settings extends AppCompatActivity {

	boolean showeda;

	// Creates the layout
	protected void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		setContentView(R.layout.settings);
		CheckBox eda = findViewById(R.id.showeda);
		CheckBox bvp = findViewById(R.id.showbvp);


		eda.setChecked(false);
		bvp.setChecked(true);

	}

	public void showeda_clicked(View v) {
		CheckBox eda = findViewById(R.id.showeda);
		CheckBox bvp = findViewById(R.id.showbvp);
		if (showeda) {
			showeda = false;
			eda.setChecked(false);
			bvp.setChecked(true);
		} else {
			showeda = true;
			eda.setChecked(true);
			bvp.setChecked(false);
		}
	}

	public void showbvp_clicked(View v){
		CheckBox eda = findViewById(R.id.showeda);
		CheckBox bvp = findViewById(R.id.showbvp);
		if(showeda){
			showeda = false;
			eda.setChecked(false);
			bvp.setChecked(true);
		}else{
			showeda = true;
			eda.setChecked(true);
			bvp.setChecked(false);
		}
	}


	public void updateGraph(View v){
		/*Intent i = new Intent();

		i.putExtra("data", showeda);
		setResult(RESULT_OK,i);
		finish();*/
	}


}