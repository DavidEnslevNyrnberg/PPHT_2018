package com.example.erik.acutestressapp;



import android.content.Intent;
import android.os.Bundle;
import android.support.v7.app.AppCompatActivity;
import android.view.View;
import android.widget.ArrayAdapter;
import android.widget.CheckBox;
import android.widget.ListView;


public class settings extends AppCompatActivity {

	// Creates the layout
	protected void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		setContentView(R.layout.settings);
		CheckBox eda = findViewById(R.id.showeda);
		CheckBox bvp = findViewById(R.id.showbvp);
		CheckBox live = findViewById(R.id.livemode);
		eda.setChecked(false);
		bvp.setChecked(true);
		live.setChecked(true);






	}

	/*public void savesettings(View v){
		Intent i = new Intent();
		String message = editText.getText().toString();
		i.putExtra("data",message);
		setResult(RESULT_OK,i);
		finish();
	}*/
}