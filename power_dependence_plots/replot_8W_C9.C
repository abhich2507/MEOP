

void replot_8W_C9() {

    // ---- 0. Open file and retrieve objects --------------------
    TFile *f = TFile::Open("8W_C9_datagrab1_four_phase_fit.root");
    if (!f || f->IsZombie()) { printf("Cannot open file!\n"); return; }

    TCanvas *cold = (TCanvas*)f->Get("c_34");
    if (!cold) { printf("Canvas c_34 not found!\n"); return; }

    // Extract graph and fits by name from the old canvas
    TGraph    *gr      = (TGraph*)cold->GetPrimitive("Graph");
    TF1       *f_ic    = (TF1*)  cold->GetPrimitive("fit_ic_34");
    TF1       *f_bu    = (TF1*)  cold->GetPrimitive("fit_bu_34");
    TF1       *f_pl    = (TF1*)  cold->GetPrimitive("fit_pl_34");
    TF1       *f_dc    = (TF1*)  cold->GetPrimitive("fit_dc_34");

    if (!gr) { printf("TGraph not found!\n"); return; }

    // ---- 1. Read fit parameters needed for the header ---------
    //  P_plateau from plateau constant fit (parameter 0)
    double P_max      = f_pl ? f_pl->GetParameter(0) : 0;
    double P_max_err  = f_pl ? f_pl->GetParError(0)  : 0;
    //  tau_b from buildup fit (parameter 2)
    double tau_b      = f_bu ? f_bu->GetParameter(2) : 0;
    double tau_b_err  = f_bu ? f_bu->GetParError(2)  : 0;
    //  tau (decay) from decay fit (parameter 1)
    double tau_r      = f_dc ? f_dc->GetParameter(1) : 0;
    double tau_r_err  = f_dc ? f_dc->GetParError(1)  : 0;

    // ---- 2. New canvas ----------------------------------------
    TCanvas *c = new TCanvas("c_clean", "8W C9 datagrab1 – Clean", 1300, 800);
    c->SetGrid();
    c->SetLeftMargin(0.10);
    c->SetRightMargin(0.03);
    c->SetTopMargin(0.17);   // extra room at top for the header text
    c->SetBottomMargin(0.12);

    // ---- 3. Draw data -----------------------------------------
    gr->SetTitle("");        // suppress ROOT title box
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.5);
    gr->SetMarkerColor(kBlack);
    gr->Draw("AP");

    // Axis cosmetics
    TAxis *xax = gr->GetXaxis();
    TAxis *yax = gr->GetYaxis();

    xax->SetTitle("Time (s)");
    yax->SetTitle("Polarization (%)");

    xax->SetTitleSize(0.055);
    yax->SetTitleSize(0.055);
    xax->SetTitleOffset(1.0);
    yax->SetTitleOffset(0.85);

    xax->SetLabelSize(0.045);
    yax->SetLabelSize(0.045);

    xax->CenterTitle(kTRUE);
    yax->CenterTitle(kTRUE);

    // ---- 4. Draw fit curves -----------------------------------
    struct FitInfo { TF1 *fn; Color_t col; };
    FitInfo fits[] = {
        {f_ic, kBlue   },
        {f_bu, kGreen+2},
        {f_pl, kOrange+7},
        {f_dc, kRed    }
    };
    for (auto &fi : fits) {
        if (!fi.fn) continue;
        fi.fn->SetLineColor(fi.col);
        fi.fn->SetLineWidth(3);
        fi.fn->Draw("SAME");
    }

    // // ---- 5. Legend --------------------------------------------
    // TLegend *leg = new TLegend(0.12, 0.62, 0.38, 0.83);
    // leg->SetTextSize(0.035);
    // leg->SetBorderSize(1);
    // leg->SetFillColor(kWhite);
    // leg->AddEntry(gr,   "Data",                  "p");
    // if (f_ic) leg->AddEntry(f_ic, "Initial constant", "l");
    // if (f_bu) leg->AddEntry(f_bu, "Buildup",          "l");
    // if (f_pl) leg->AddEntry(f_pl, "Plateau",          "l");
    // if (f_dc) leg->AddEntry(f_dc, "Decay",            "l");
    // leg->Draw();

    // ---- 6. Header text (NDC coordinates, top of pad) ---------
    // Three lines stacked vertically, all centred at x0
    double x0    = 0.50;   // horizontal centre of canvas
    double y_top = 0.975;  // first line
    double dy    = 0.052;  // line spacing
    double sz    = 0.042;  // text size (slightly bigger to be readable)

    auto makeLabel = [&](double x, double y, const char *txt) {
        TLatex *l = new TLatex(x, y, txt);
        l->SetNDC();
        l->SetTextSize(sz);
        l->SetTextAlign(22);  // centre-centre
        l->SetTextFont(42);
        l->Draw();
        return l;
    };

    // Line 1: P_max
    TString sP = TString::Format(
        "P_{max} = %.2f #pm %.2f %%", P_max, P_max_err);
    makeLabel(x0, y_top, sP);

    // Line 2: T_pump (tau_b)
    TString sTb = TString::Format(
        "T_{pump} (#tau_{b}) = %.2f #pm %.2f s", tau_b, tau_b_err);
    makeLabel(x0, y_top - dy, sTb);

    // Line 3: T_relax (tau_decay)
    TString sTr = TString::Format(
        "T_{relax} (#tau_{decay}) = %.2f #pm %.2f s", tau_r, tau_r_err);
    makeLabel(x0, y_top - 2*dy, sTr);

    // ---- 7. Save ----------------------------------------------
    c->Update();
    c->SaveAs("8W_C9_datagrab1_clean.png");
    c->SaveAs("8W_C9_datagrab1_clean.pdf");
    printf("\nSaved: 8W_C9_datagrab1_clean.png / .pdf\n");
    printf("  P_max   = %.3f +/- %.3f %%\n",  P_max,  P_max_err);
    printf("  tau_b   = %.3f +/- %.3f s\n",   tau_b,  tau_b_err);
    printf("  tau_rel = %.3f +/- %.3f s\n",   tau_r,  tau_r_err);
}
